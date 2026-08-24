package webserver

import (
	"os"
	"path/filepath"
	"testing"
)

// LoadUsers must not write when the database file is absent.
//
// init() calls LoadUsers, and the path it uses is relative to the PROCESS working
// directory. A write there means merely importing this package creates a users database
// wherever it happens to be running — which left an untracked pkg/webserver/users/users.json
// behind on every `go test ./pkg/webserver/`, and would shadow the real database for any
// tool started from outside the DRIFT folder.
func TestLoadUsersDoesNotWriteWhenFileMissing(t *testing.T) {
	orig, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	if err := os.Chdir(dir); err != nil {
		t.Fatal(err)
	}
	// Restore before t.TempDir's cleanup runs, or the directory cannot be removed.
	defer func() {
		if err := os.Chdir(orig); err != nil {
			t.Fatalf("restoring working directory: %v", err)
		}
	}()

	db, err := LoadUsers()
	if err != nil {
		t.Fatalf("LoadUsers with no file present: %v", err)
	}
	if len(db.Users) == 0 {
		t.Error("expected the in-memory placeholder user to be returned")
	}
	if _, err := os.Stat(filepath.Join(dir, "users")); !os.IsNotExist(err) {
		t.Error("LoadUsers created a users/ directory; reading the database must not have " +
			"a write side effect, because init() calls it from whatever cwd the process has")
	}
}

// The placeholder is only in memory, so a real registration is what creates the file.
func TestAddUserPersists(t *testing.T) {
	orig, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	if err := os.Chdir(dir); err != nil {
		t.Fatal(err)
	}
	defer func() {
		if err := os.Chdir(orig); err != nil {
			t.Fatalf("restoring working directory: %v", err)
		}
	}()

	db, err := LoadUsers()
	if err != nil {
		t.Fatalf("LoadUsers: %v", err)
	}
	if err := db.AddUser("rob", ""); err != nil {
		t.Fatalf("AddUser: %v", err)
	}
	if _, err := os.Stat(filepath.Join(dir, "users", "users.json")); err != nil {
		t.Errorf("registration did not persist the database: %v", err)
	}
}
