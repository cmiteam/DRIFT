package checkpoint

import (
	"testing"
)

func TestStopRequestLifecycle(t *testing.T) {
	dir := t.TempDir()

	if _, ok := StopRequested(dir); ok {
		t.Fatal("no stop should be requested initially")
	}

	if err := RequestStop(dir, "baseline"); err != nil {
		t.Fatalf("RequestStop: %v", err)
	}
	name, ok := StopRequested(dir)
	if !ok || name != "baseline" {
		t.Errorf("StopRequested = (%q, %v), want (baseline, true)", name, ok)
	}

	ClearStop(dir)
	if _, ok := StopRequested(dir); ok {
		t.Error("stop request should be cleared")
	}
}

func TestRequestStopEmptyName(t *testing.T) {
	dir := t.TempDir()
	if err := RequestStop(dir, ""); err != nil {
		t.Fatalf("RequestStop: %v", err)
	}
	name, ok := StopRequested(dir)
	if !ok || name != "" {
		t.Errorf("StopRequested = (%q, %v), want (\"\", true)", name, ok)
	}
}
