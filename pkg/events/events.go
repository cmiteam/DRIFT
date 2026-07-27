package events

import (
	"drift/pkg/core"
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// Scriptable environmental events (roadmap §2).
//
// An event is a scheduled perturbation of the population — a famine that raises
// mortality for a span of years, a one-time migration pulse between demes, and so
// on. Events generalize the old single-purpose events package (which only held the
// seeding modules): each event TYPE self-registers under a string key (mirroring
// the pkg/modules phase registry) and knows how to (a) parse its own key=value
// fields and (b) apply its effect for a given year.
//
// A model opts in via the `environmental_events` string param. LoadSchedule parses
// it into a Schedule (implementing core.EventScheduler) which the run loop calls
// once per year. When the param is empty (the default) no scheduler is attached,
// so the engine takes its byte-identical existing path — the strict no-op the §6h
// neutral-validation harness depends on. Because the schedule is a pure function of
// the param string (no serialized state), checkpoints need no schema bump: it is
// rebuilt from the params on resume, exactly like the demographic scheduler.

// Event is one scheduled environmental perturbation. Apply is called once per
// simulated year with the current year; the event decides for itself whether it
// acts this year (a windowed famine checks its span; a one-time migration pulse
// checks its year). An event MUST consume RNG only in sorted-id order, and none at
// all in years it does not act, so a fixed rng_seed stays byte-reproducible and
// years before the first event match a no-event baseline exactly.
type Event interface {
	Apply(model *core.Model, pop *core.Pop, year int)
}

// Factory builds an event of a given type from its parsed key=value fields,
// returning a clear error for missing/invalid fields (surfaced at model load).
type Factory func(fields map[string]string) (Event, error)

var eventRegistry = map[string]Factory{}

// RegisterEvent registers an event type under a string key. Call it from an init()
// in the file that implements the type (mirrors pkg/modules' phase registries).
func RegisterEvent(name string, f Factory) { eventRegistry[name] = f }

// AvailableEvents returns the registered event-type names, sorted — for error
// messages and a future GUI/list endpoint.
func AvailableEvents() []string {
	names := make([]string, 0, len(eventRegistry))
	for n := range eventRegistry {
		names = append(names, n)
	}
	sort.Strings(names)
	return names
}

// Schedule is a parsed, ordered list of events. It implements core.EventScheduler.
type Schedule struct {
	events []Event
}

// Len reports how many events the schedule holds (for load-time logging).
func (s *Schedule) Len() int { return len(s.events) }

// Apply advances every scheduled event for the current year. It first resets the
// per-year event modifiers on the model to neutral (so a famine that ended last
// year stops raising mortality), then applies each event in the order it was
// declared. Events run in that fixed spec order and each consumes RNG only in
// sorted-id order, so the whole run stays byte-reproducible under a fixed seed.
func (s *Schedule) Apply(model *core.Model, pop *core.Pop) {
	year := model.FreeParameters["year"]
	// Neutral defaults for the transient per-year modifiers, rewritten by any
	// active event below. Keeping this here (not in each event) guarantees a
	// modifier never lingers past the event that set it.
	model.EventMortalityFactor = 1.0
	for _, e := range s.events {
		e.Apply(model, pop, year)
	}
}

// LoadSchedule parses the `environmental_events` param into a Schedule, or returns
// (nil, nil) when the param is empty (the default) so the run takes the strict
// no-op path. The spec mirrors mutation_classes: a ';'-separated list of entries,
// each
//
//	<type> key=value key=value ...
//
// whitespace-tokenized (never commas, so the value needs no CSV quoting). The
// first token selects the registered event type; the remaining fields are passed
// to its Factory. An unknown type or a malformed/invalid field is a hard error
// (fail fast at model load, matching the ValidateStyles philosophy) rather than a
// silently skipped event.
//
// Example:
//
//	famine start=200 end=205 mortality=2.0; migration_pulse year=300 from=0 to=1 fraction=0.1
func LoadSchedule(model *core.Model) (*Schedule, error) {
	spec := strings.TrimSpace(model.StringParam("environmental_events", ""))
	if spec == "" {
		return nil, nil
	}
	var events []Event
	for _, entry := range strings.Split(spec, ";") {
		entry = strings.TrimSpace(entry)
		if entry == "" {
			continue
		}
		ev, err := parseEvent(entry)
		if err != nil {
			return nil, err
		}
		events = append(events, ev)
	}
	if len(events) == 0 {
		return nil, nil
	}
	return &Schedule{events: events}, nil
}

// parseEvent turns one spec entry into an Event via its registered factory.
func parseEvent(entry string) (Event, error) {
	tokens := strings.Fields(entry)
	if len(tokens) == 0 {
		return nil, fmt.Errorf("events: empty event entry")
	}
	name := tokens[0]
	factory, ok := eventRegistry[name]
	if !ok {
		return nil, fmt.Errorf("events: unknown event type %q; available: %v", name, AvailableEvents())
	}
	fields := make(map[string]string, len(tokens)-1)
	for _, t := range tokens[1:] {
		kv := strings.SplitN(t, "=", 2)
		if len(kv) != 2 {
			return nil, fmt.Errorf("events: %s: malformed field %q (want key=value)", name, t)
		}
		fields[strings.ToLower(strings.TrimSpace(kv[0]))] = strings.TrimSpace(kv[1])
	}
	return factory(fields)
}

// Field-parsing helpers shared by the event factories.

func reqInt(fields map[string]string, key string) (int, error) {
	v, ok := fields[key]
	if !ok {
		return 0, fmt.Errorf("missing required field %q", key)
	}
	n, err := strconv.Atoi(v)
	if err != nil {
		return 0, fmt.Errorf("field %q: %q is not an integer", key, v)
	}
	return n, nil
}

func optInt(fields map[string]string, key string, def int) (int, error) {
	if _, ok := fields[key]; !ok {
		return def, nil
	}
	return reqInt(fields, key)
}

func reqFloat(fields map[string]string, key string) (float64, error) {
	v, ok := fields[key]
	if !ok {
		return 0, fmt.Errorf("missing required field %q", key)
	}
	f, err := strconv.ParseFloat(v, 64)
	if err != nil {
		return 0, fmt.Errorf("field %q: %q is not a number", key, v)
	}
	return f, nil
}
