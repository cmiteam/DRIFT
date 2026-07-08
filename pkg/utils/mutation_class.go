package utils

import (
	"math"
	"strconv"
	"strings"

	"drift/pkg/core"
)

// MutationClass describes one class in the mutation spectrum (roadmap §1). Real
// genomes accumulate several distinct kinds of mutation — point substitutions,
// small indels, copy-number variants, large deletions — each with its own
// per-reproduction rate and its own distribution of fitness effects. DRIFT
// historically drew a single uniform mu with one Weibull DFE; a MutationClass
// lets each class carry independent draws so the spectrum is configurable and,
// once recorded on each Mutation, downstream statistics can partition by class.
//
// The effect machinery is deliberately the existing Weibull DFE (shape/scale/adj
// + f_neutral + f_beneficial) applied per class, so no new distribution code is
// introduced — a class differs from the legacy behavior only in its parameters.
type MutationClass struct {
	Name        string  // human label (point / indel / CNV / large_deletion / ...)
	Rate        float64 // Poisson mean of new mutations of this class per reproduction
	FNeutral    float64 // fraction of this class that is neutral (Effect stays 0)
	FBeneficial float64 // among non-neutral, fraction beneficial (else deleterious)
	Shape       float64 // Weibull shape for the |effect| magnitude
	Scale       float64 // Weibull scale
	WeibullAdj  float64 // divisor applied to the drawn magnitude
	Dominance   int      // h*100 stamped on new mutations of this class
	Size        int      // target size in genome bits recorded on each mutation
}

// dominancePct reads the global dominance coefficient h (fitness_dominance) as an
// int percentage. Default 50 (additive) when the param is absent — a bare map
// lookup would yield h=0 (fully recessive), which is not the intended default.
func dominancePct(model *core.Model) int {
	pct := 50
	if h, ok := model.Parameters["fitness_dominance"]; ok {
		pct = int(math.Round(h * 100))
	}
	return pct
}

// defaultMutationClass is the legacy single class: a "point" class whose rate and
// effect distribution are the model's global mu / Weibull / dominance settings.
// When mutation_classes is unset this is the ONLY class, and GenerateNewMutations
// then draws exactly the same RNG sequence as before the spectrum existed (one
// Poisson(mu), then the same per-mutation draws) — the strict no-op the §6h
// neutral-validation harness depends on.
func defaultMutationClass(model *core.Model) MutationClass {
	return MutationClass{
		Name:        "point",
		Rate:        model.Parameters["mu"],
		FNeutral:    model.Parameters["f_neutral"],
		FBeneficial: model.Parameters["f_beneficial"],
		Shape:       model.Parameters["shape"],
		Scale:       model.Parameters["scale"],
		WeibullAdj:  model.Parameters["Weibull_adj"],
		Dominance:   dominancePct(model),
		Size:        1,
	}
}

// MutationClasses returns the configured mutation-class spectrum.
//
// With the mutation_classes string param empty (the default) it returns a single
// point class (defaultMutationClass) — byte-identical legacy behavior. A non-empty
// spec is a ';'-separated list of classes; each class is whitespace-tokenized as
//
//	<name> key=value key=value ...
//
// where the first token is the class name and the remaining key=value fields
// override the inherited defaults. Recognized keys: rate, fneutral, fbeneficial,
// shape, scale, adj (weibull_adj), dom (h in 0..1), size. Any field omitted for a
// class inherits the model's global value, so e.g. `indel rate=1 scale=0.2` reuses
// the global Weibull shape / f_neutral / f_beneficial. Fields use whitespace and
// ';' separators (never commas) so the value needs no CSV quoting.
//
// Example (point substitutions + rarer, larger-effect indels and large deletions):
//
//	point rate=8; indel rate=1 scale=0.15 size=10; large_deletion rate=0.1 fneutral=0 fbeneficial=0 scale=0.4 size=1000
func MutationClasses(model *core.Model) []MutationClass {
	spec := strings.TrimSpace(model.StringParam("mutation_classes", ""))
	base := defaultMutationClass(model)
	if spec == "" {
		return []MutationClass{base}
	}

	var classes []MutationClass
	for _, entry := range strings.Split(spec, ";") {
		entry = strings.TrimSpace(entry)
		if entry == "" {
			continue
		}
		classes = append(classes, parseMutationClass(entry, base))
	}
	if len(classes) == 0 {
		// Spec was all separators/whitespace — fall back to the legacy no-op.
		return []MutationClass{base}
	}
	return classes
}

// parseMutationClass parses one class entry, starting from base (the model
// globals) and overriding only the fields the entry names.
func parseMutationClass(entry string, base MutationClass) MutationClass {
	fields := strings.Fields(entry)
	c := base
	if len(fields) == 0 {
		return c
	}
	c.Name = fields[0]
	for _, f := range fields[1:] {
		kv := strings.SplitN(f, "=", 2)
		if len(kv) != 2 {
			continue
		}
		key := strings.ToLower(strings.TrimSpace(kv[0]))
		val := strings.TrimSpace(kv[1])
		switch key {
		case "rate":
			c.Rate = parseFloatDefault(val, c.Rate)
		case "fneutral", "f_neutral":
			c.FNeutral = parseFloatDefault(val, c.FNeutral)
		case "fbeneficial", "f_beneficial":
			c.FBeneficial = parseFloatDefault(val, c.FBeneficial)
		case "shape":
			c.Shape = parseFloatDefault(val, c.Shape)
		case "scale":
			c.Scale = parseFloatDefault(val, c.Scale)
		case "adj", "weibull_adj":
			c.WeibullAdj = parseFloatDefault(val, c.WeibullAdj)
		case "dom", "dominance":
			h := parseFloatDefault(val, float64(c.Dominance)/100.0)
			c.Dominance = int(math.Round(h * 100))
		case "size":
			c.Size = int(math.Round(parseFloatDefault(val, float64(c.Size))))
		}
	}
	return c
}

func parseFloatDefault(s string, def float64) float64 {
	if v, err := strconv.ParseFloat(s, 64); err == nil {
		return v
	}
	return def
}
