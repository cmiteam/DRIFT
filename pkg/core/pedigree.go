package core

// Diploid pedigree retention — TMR4A.md work item W1 (`track_pedigree`).
//
// WHY THIS EXISTS. DRIFT already records two ancestry chains: MaleDB (child -> father,
// written only when the child is male) and FemaleDB (child -> mother, written only when
// the child is female). Those are exactly what FindYAdam / FindMtEve need — the Y-chain
// and the mt-chain — and they are NOT a pedigree: a son's mother and a daughter's father
// are recorded nowhere, so no autosomal lineage can be traced. The TMR-K-A walker (W5)
// walks (individual, strand) pairs backwards through ordinary autosomal inheritance, so
// it needs BOTH parents of EVERY child, and it needs them to survive the parent's death.
//
// THE RETENTION RULE. A node must be kept for as long as any sampled descendant might
// still walk through it, and dropped the moment it cannot be reached. That is a refcount:
//
//	Refs = (1 if this individual is still alive) + (number of RETAINED child nodes)
//
// A node is released when Refs hits zero, and releasing it decrements its parents, which
// may release them in turn — the cascade prunes whole extinct branches in one pass. This
// generalises the existing MaleDB/FemaleDB prune rule in death.go ("drop a dead node only
// when it left no same-sex offspring"), which is the same discipline specialised to a
// single-sex chain. It is also the same discipline as the §6f mutation-pool refcount, and
// it has the same failure mode: an asymmetric increment/decrement leaks nodes forever or
// frees a node someone still points at. Every increment here has exactly one matching
// decrement, and the pairing is asserted by TestPedigreeRefcountBalance.
//
// WHY BIRTHYEAR AND SEX ARE COPIED IN. IndData[ind] is deleted unconditionally at
// death.go:325, so after a parent dies the only surviving record of when it was born is
// whatever the pedigree kept. The walker reports coalescence in simulated YEARS, so the
// birth year has to live here. Sex is carried because it costs one int and lets the
// walker report a genealogy without a second lookup that would fail for the dead.
//
// COST AND DETERMINISM. This path draws NO RNG and touches no existing statistic, so a
// track_pedigree run is byte-identical to a track_pedigree-off run (the §0 gate). It is
// opt-in: with track_pedigree unset, Pop.Pedigree stays nil and every method below is a
// no-op returning immediately. Memory is the real cost and the item most likely to force
// a design change at mainstream scale (TMR4A.md §5) — one node is ~6 words plus map
// overhead, retained for everyone who left surviving descendants, so measure it early.

// PedNode is one retained individual in the diploid pedigree. Exported (as are its
// fields) so the whole structure survives gob checkpointing with no custom codec.
type PedNode struct {
	Dad       int  // father's ID, or PedFounder for a founder / unrecorded parent
	Mom       int  // mother's ID, or PedFounder
	BirthYear int  // simulated year of birth (negative for pre-simulation founders)
	Sex       int  // 0 = male, 1 = female
	Refs      int  // retained children + 1 while alive; released at 0
	Alive     bool // false once RIP has run, so the alive-contribution is removed exactly once
}

// PedFounder marks a parent slot with no recorded ancestor — a created/initial founder,
// or an individual born before track_pedigree was switched on. The walker stops here and
// reports `censored_at_founding` rather than inventing a coalescence (TMR4A.md §0: the
// backward walk must terminate at created founders and SAY SO, never force a root).
const PedFounder = -1

// PedigreeOn reports whether pedigree retention is active for this population.
func (p *Pop) PedigreeOn() bool { return p != nil && p.Pedigree != nil }

// PedEnsure creates a founder-style node for id if none exists yet, and returns it.
// Used for the initial population (seeded once after setup) and defensively for any
// parent that turns out to be untracked — an individual created outside birth.go, or a
// checkpoint resumed with track_pedigree newly enabled. A node created here has no
// recorded parents, so the walk terminates at it, which is the correct and honest answer
// for someone whose ancestry was never recorded.
//
// The node starts Alive with Refs = 1; RIP removes that contribution exactly once.
func (p *Pop) PedEnsure(id, birthYear, sex int) *PedNode {
	if p.Pedigree == nil {
		p.Pedigree = make(map[int]*PedNode)
	}
	if n, ok := p.Pedigree[id]; ok {
		return n
	}
	n := &PedNode{Dad: PedFounder, Mom: PedFounder, BirthYear: birthYear, Sex: sex, Refs: 1, Alive: true}
	p.Pedigree[id] = n
	return n
}

// PedRecordBirth records a child with BOTH parents — the whole point of W1 — and takes a
// reference on each parent that exists. Parents are expected to have been ensured by the
// caller (birth.go ensures them from IndData, which is still live for a parent); a parent
// with no node is simply not referenced, and the child's slot still names it, so a walk
// that reaches a missing parent terminates as censored rather than panicking.
func (p *Pop) PedRecordBirth(child, dad, mom, birthYear, sex int) {
	if p.Pedigree == nil {
		p.Pedigree = make(map[int]*PedNode)
	}
	p.Pedigree[child] = &PedNode{Dad: dad, Mom: mom, BirthYear: birthYear, Sex: sex, Refs: 1, Alive: true}
	if n, ok := p.Pedigree[dad]; ok {
		n.Refs++
	}
	if n, ok := p.Pedigree[mom]; ok {
		n.Refs++
	}
	// Self-mating is impossible in DRIFT (dad is male, mom is female), so the two
	// increments above can never both land on the same node.
}

// PedRecordDeath removes the alive-contribution of id and prunes the branch if that was
// the last reference. Called once per individual, from RIP, BEFORE IndData is deleted.
// Idempotent in the sense that a second call on an already-dead node does nothing: the
// Alive flag guards the decrement, so a double-RIP cannot free a node someone points at.
func (p *Pop) PedRecordDeath(id int) {
	if p.Pedigree == nil {
		return
	}
	n, ok := p.Pedigree[id]
	if !ok || !n.Alive {
		return
	}
	n.Alive = false
	n.Refs--
	p.pedRelease(id)
}

// pedRelease drops id if it is unreferenced and cascades to its parents. Iterative rather
// than recursive on purpose: a lineage can be tens of thousands of generations deep in a
// long run, and a recursive cascade would blow the stack on the one run that matters.
func (p *Pop) pedRelease(id int) {
	stack := []int{id}
	for len(stack) > 0 {
		cur := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		n, ok := p.Pedigree[cur]
		if !ok || n.Refs > 0 || n.Alive {
			continue
		}
		dad, mom := n.Dad, n.Mom
		delete(p.Pedigree, cur)
		// Strand provenance (W2) is indexed BY the pedigree and is meaningless without
		// it, so it is pruned by the same cascade rather than by a second, separately
		// maintained refcount that could drift out of step with this one.
		delete(p.ARGDB, cur)

		for _, parent := range [2]int{dad, mom} {
			if parent == PedFounder {
				continue
			}
			pn, ok := p.Pedigree[parent]
			if !ok {
				continue
			}
			pn.Refs--
			if pn.Refs <= 0 && !pn.Alive {
				stack = append(stack, parent)
			}
		}
	}
}

// PedAncestryPath returns the chain of nodes from id back to a founder along a single
// named parent slot, newest first. Diagnostic helper for tests and reporting; the W5
// walker does its own multi-lineage traversal and does not use this.
func (p *Pop) PedAncestryPath(id int, maternal bool) []int {
	if p.Pedigree == nil {
		return nil
	}
	var path []int
	for cur := id; cur != PedFounder; {
		n, ok := p.Pedigree[cur]
		if !ok {
			break
		}
		path = append(path, cur)
		if maternal {
			cur = n.Mom
		} else {
			cur = n.Dad
		}
	}
	return path
}
