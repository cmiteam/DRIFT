package webserver

import (
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

// sanitizeOutputDir guards a filesystem WRITE target supplied by a web form, so the
// cases that matter are the ones that try to leave the DRIFT folder. An accepted value
// must always come back relative and inside the tree; anything else must be an error
// rather than a surprising write.

const fallback = "users/rob/models/m/results"

func TestSanitizeOutputDirEmptyKeepsDefault(t *testing.T) {
	for _, in := range []string{"", "   ", ".", "./"} {
		got, err := sanitizeOutputDir(in, fallback)
		if err != nil {
			t.Errorf("%q: unexpected error %v", in, err)
			continue
		}
		if got != fallback {
			t.Errorf("%q: got %q, want the per-model default %q", in, got, fallback)
		}
	}
}

func TestSanitizeOutputDirAcceptsRelativePaths(t *testing.T) {
	root, err := filepath.Abs(".")
	if err != nil {
		t.Fatal(err)
	}
	for _, in := range []string{
		"results/argweaver/createdvar",
		`results\argweaver\createdvar`, // Windows separators must work too
		"results/../results/ok",        // interior .. that does not escape
		"results",
	} {
		got, err := sanitizeOutputDir(in, fallback)
		if err != nil {
			t.Errorf("%q: unexpected error %v", in, err)
			continue
		}
		if filepath.IsAbs(got) {
			t.Errorf("%q: returned an absolute path %q", in, got)
		}
		abs, _ := filepath.Abs(filepath.Join(root, filepath.FromSlash(got)))
		rel, err := filepath.Rel(root, abs)
		if err != nil || rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
			t.Errorf("%q: resolved to %q, which is outside the DRIFT folder", in, got)
		}
	}
}

func TestSanitizeOutputDirRejectsEscapes(t *testing.T) {
	for _, in := range []string{
		"..",
		"../outside",
		"results/../../outside",
		`..\outside`,
		"results/../../../../../../etc",
	} {
		got, err := sanitizeOutputDir(in, fallback)
		if err == nil {
			t.Errorf("%q was accepted as %q; a path that climbs out of the DRIFT folder "+
				"must be refused, not written to", in, got)
		}
	}
}

func TestSanitizeOutputDirRejectsAbsolute(t *testing.T) {
	cases := []string{filepath.Join(os.TempDir(), "drift-escape")}
	if runtime.GOOS == "windows" {
		// "C:foo" is drive-RELATIVE: filepath.IsAbs reports false for it, so the volume
		// check is what catches it. UNC paths are the other Windows-only escape.
		cases = append(cases, `C:foo`, `C:\Windows\Temp`, `\\server\share\x`)
	}
	for _, in := range cases {
		got, err := sanitizeOutputDir(in, fallback)
		if err == nil {
			t.Errorf("%q was accepted as %q; absolute and volume-qualified paths must be "+
				"refused", in, got)
		}
	}
}
