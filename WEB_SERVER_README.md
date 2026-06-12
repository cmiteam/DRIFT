# Drift Web Server - Implementation Summary

## Overview
The Drift web server has been successfully implemented to prepare for the HTML GUI. The server provides REST API endpoints for user authentication, simulation management, and result access.

## Running the Web Server

### Start the Server
```bash
./drift.exe --web
```

The server will start on `http://localhost:8080`

### CLI Mode (Still Works!)
```bash
./drift.exe                    # Run simulation with defaults
./drift.exe --results my_dir   # Custom results directory
```

## Architecture

### Package Structure
```
pkg/webserver/
├── server.go     # Main HTTP server setup
├── auth.go       # User authentication
├── queue.go      # Simulation queue management
└── handlers.go   # API endpoint handlers
```

### Directory Structure
```
users/
├── users.json           # User authentication database
├── alice/
│   └── Default/
│       ├── progress.json       # Real-time progress updates
│       ├── Default_results.csv
│       ├── genetic.gif
│       └── genealogical.gif
└── bob/
    └── Default/
        └── ...

web/static/
└── index.html          # Placeholder GUI (to be expanded)
```

## API Endpoints

### Authentication
- **POST** `/api/login`
  - Body: `{"username": "alice", "password": ""}`
  - Response: `{"success": true, "username": "alice"}`

- **GET** `/api/users`
  - Returns list of all usernames

### Simulation Management
- **POST** `/api/simulation/start`
  - Body: `{"username": "alice", "model_name": "Default"}`
  - Queues a new simulation
  - Response: Job object with ID and status

- **GET** `/api/simulation/status?username=alice`
  - Returns all jobs for a user

- **POST** `/api/simulation/cancel`
  - Body: `{"job_id": "alice_Default_1234567890"}`
  - Cancels queued or running simulation

### Progress & Results
- **GET** `/api/progress/alice/Default`
  - Returns the current progress.json file

- **GET** `/api/results/alice/Default`
  - Lists all result files

- **GET** `/api/results/alice/Default/genetic.gif`
  - Downloads a specific result file

### Admin (Future: Add authentication)
- **GET** `/api/admin/queue`
  - View full simulation queue

- **GET** `/api/admin/users`
  - View all users

- **POST** `/api/admin/users`
  - Create new user

## Features Implemented

### ✅ User Management
- User authentication via users.json
- Multiple users supported
- Password optional (empty string = no password)

### ✅ Simulation Queue
- First-come-first-served queue
- Only one simulation runs at a time (local deployment)
- Users can cancel their own simulations
- Real-time status tracking

### ✅ Process Management
- Web server spawns `drift.exe` as subprocess
- Command: `drift.exe --output-dir users/{username}/{modelname}`
- Monitors process completion
- Handles errors and cancellations

### ✅ Progress Tracking
- `progress.json` updated at each save_interval
- Shows current generation, population size, etc.
- Frontend can poll for updates

### ✅ Results Access
- List all result files via API
- Download CSVs, GIFs, PNGs
- Security: Path traversal protection

## How It Works

### User Workflow
1. User opens browser to `http://localhost:8080`
2. Logs in with username (password optional)
3. Configures simulation parameters (GUI to be built)
4. Clicks "Run Simulation"
5. Job added to queue
6. When job starts, drift.exe spawns
7. Frontend polls `/api/progress/{username}/{model}` for updates
8. When complete, view results

### Behind the Scenes
```
Browser --> HTTP Request --> Web Server --> Queue --> drift.exe subprocess
                                  ↓                         ↓
                            users.json              users/alice/Default/
                                                    ├── progress.json
                                                    ├── results.csv
                                                    └── *.gif
```

## Testing

### Test Login API
```bash
curl -X POST http://localhost:8080/api/login \
  -H "Content-Type: application/json" \
  -d '{"username":"alice","password":""}'
```

### Test Start Simulation
```bash
curl -X POST http://localhost:8080/api/simulation/start \
  -H "Content-Type: application/json" \
  -d '{"username":"alice","model_name":"Default"}'
```

### Test Progress
```bash
curl http://localhost:8080/api/progress/alice/Default
```

## Next Steps

### Immediate (Before HTML GUI)
- ✅ Web server implementation
- ⏳ Create full HTML/CSS/JavaScript frontend
- ⏳ Build 4-tab interface (Setup, Progress, Graphs, Maps)

### HTML GUI Tabs
1. **Setup Tab**
   - Form for all parameters from parameter_defaults.csv
   - Load/save configurations
   - Start simulation button

2. **Progress Tab**
   - Real-time progress updates (polling /api/progress)
   - Current generation, population size
   - Queue position if queued
   - Cancel button

3. **Graphs Tab**
   - Display genetic.gif and genealogical.gif
   - Auto-refresh as they update
   - Download buttons

4. **Maps Tab**
   - Display map animations (if enabled)
   - Conditional on track_map parameter

## Security Notes

### Current State (Local Development)
- No password encryption (TODO: bcrypt)
- No session management
- No admin authentication
- Simple path traversal protection

### For Production/Cloud
- Implement proper bcrypt password hashing
- Add session tokens/JWT
- Rate limiting
- Admin authentication
- HTTPS/TLS
- CORS configuration

## Configuration

### Environment Variables
- `PORT` - Server port (default: 8080)

### Customization
Edit `pkg/webserver/server.go` to change:
- Port number
- Static file directory
- Drift executable path

## Troubleshooting

### "drift.exe not found"
- Ensure drift.exe is in the same directory
- Or build: `go build -o drift.exe drift.go`

### "web/static not found"
- Create directory: `mkdir -p web/static`
- Add index.html

### Port already in use
- Change port: `export PORT=3000`
- Or kill process on port 8080

## Development

### Add New Endpoint
1. Create handler function in `pkg/webserver/handlers.go`
2. Register route in `pkg/webserver/server.go`
3. Document in this README

### Modify Queue Behavior
- Edit `pkg/webserver/queue.go`
- Change `processor()` for parallel execution
- Adjust priority logic

## Compatibility

- Works on Windows, macOS, Linux
- Pure Go HTTP server (no dependencies)
- Compatible with existing CLI mode
- Backward compatible with all drift.go features
