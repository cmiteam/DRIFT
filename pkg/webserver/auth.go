package webserver

import (
	"encoding/json"
	"fmt"
	"os"
	"sync"
	"time"
)

// User represents a user in the system
type User struct {
	Username     string    `json:"username"`
	PasswordHash string    `json:"password_hash"`
	Created      time.Time `json:"created"`
	Models       []string  `json:"models"`
}

// UsersDatabase represents the users.json file structure
type UsersDatabase struct {
	Users []User `json:"users"`
	mu    sync.RWMutex
}

var usersDB *UsersDatabase

// LoadUsers loads users from users.json
func LoadUsers() (*UsersDatabase, error) {
	db := &UsersDatabase{
		Users: []User{},
	}

	data, err := os.ReadFile("users/users.json")
	if err != nil {
		// If file doesn't exist, create default
		if os.IsNotExist(err) {
			db.Users = []User{
				{
					Username:     "alice",
					PasswordHash: "",
					Created:      time.Now(),
					Models:       []string{},
				},
			}
			return db, db.Save()
		}
		return nil, fmt.Errorf("failed to read users.json: %v", err)
	}

	err = json.Unmarshal(data, db)
	if err != nil {
		return nil, fmt.Errorf("failed to parse users.json: %v", err)
	}

	return db, nil
}

// Save writes the users database to disk
func (db *UsersDatabase) Save() error {
	db.mu.RLock()
	defer db.mu.RUnlock()

	data, err := json.MarshalIndent(db, "", "  ")
	if err != nil {
		return fmt.Errorf("failed to marshal users: %v", err)
	}

	err = os.WriteFile("users/users.json", data, 0644)
	if err != nil {
		return fmt.Errorf("failed to write users.json: %v", err)
	}

	return nil
}

// GetUser returns a user by username
func (db *UsersDatabase) GetUser(username string) (*User, error) {
	db.mu.RLock()
	defer db.mu.RUnlock()

	for i := range db.Users {
		if db.Users[i].Username == username {
			return &db.Users[i], nil
		}
	}

	return nil, fmt.Errorf("user not found: %s", username)
}

// AddUser adds a new user to the database
func (db *UsersDatabase) AddUser(username, passwordHash string) error {
	db.mu.Lock()
	defer db.mu.Unlock()

	// Check if user already exists
	for _, user := range db.Users {
		if user.Username == username {
			return fmt.Errorf("user already exists: %s", username)
		}
	}

	newUser := User{
		Username:     username,
		PasswordHash: passwordHash,
		Created:      time.Now(),
		Models:       []string{},
	}

	db.Users = append(db.Users, newUser)
	return db.Save()
}

// ValidateUser checks if username/password combination is valid
func (db *UsersDatabase) ValidateUser(username, password string) bool {
	user, err := db.GetUser(username)
	if err != nil {
		return false
	}

	// If no password hash is set, allow login without password
	if user.PasswordHash == "" {
		return true
	}

	// TODO: Implement proper password hashing (bcrypt)
	// For now, simple comparison (NOT SECURE - placeholder only)
	return user.PasswordHash == password
}

// Initialize users database on package init
func init() {
	var err error
	usersDB, err = LoadUsers()
	if err != nil {
		fmt.Printf("Warning: Failed to load users database: %v\n", err)
		usersDB = &UsersDatabase{Users: []User{}}
	}
}
