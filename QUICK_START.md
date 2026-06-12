# Quick Start Guide - New Model System

## What Changed

The system now automatically creates 5 base models for each user:
- **Default** - Standard population (n=1000, t=100)
- **QuickTest** - Fast testing (n=100, t=500)
- **Flood** - Biblical flood scenario (n=6)
- **Creation** - Adam & Eve (n=2)
- **OoA** - Out of Africa (n=10000)

## For Web GUI Users

### 1. Delete Old Users Directory
```bash
rm -rf users/
```

### 2. Restart Web Server
```bash
go run drift.go --web
# OR if you built the executable:
./drift.exe --web
```

### 3. Register/Create User
- Go to http://localhost:8080
- Register a new user (e.g., "Rob")
- **All 5 base models are automatically copied to your user directory**

### 4. Select Model
- In the web GUI, you should now see a dropdown with all 5 models
- Select the model you want to run
- Click "Start Simulation"

### 5. Save/Load/Reset (NEW!)
The web GUI now supports (via API):
- **Save**: Save your current parameter changes
- **Load**: Load a specific model's parameters
- **Reset**: Reset model back to base defaults

## Directory Structure

After creating user "Rob", you'll have:
```
users/Rob/models/
├── Default/
│   ├── metadata.json
│   ├── config/parameters.csv
│   └── results/
├── QuickTest/
├── Flood/
├── Creation/
└── OoA/
```

## Testing from Command Line

```bash
# Initialize a new user
go run test_user_init.go

# Run with a specific model
./drift.exe -username=Rob -model=Default

# Run quick test
./drift.exe -username=Rob -model=QuickTest
```

## API Endpoints (for frontend integration)

```bash
# List user's models
curl "http://localhost:8080/api/models/list?username=Rob"

# Load model parameters
curl "http://localhost:8080/api/models/load?username=Rob&model_name=Default"

# Save model parameters
curl -X POST http://localhost:8080/api/models/save \
  -H "Content-Type: application/json" \
  -d '{
    "username": "Rob",
    "model_name": "Default",
    "parameters": {
      "end_year": 1000,
      "start_pop_size": 500
    }
  }'

# Reset model to defaults
curl -X POST http://localhost:8080/api/models/reset \
  -H "Content-Type: application/json" \
  -d '{
    "username": "Rob",
    "model_name": "Default"
  }'
```

## Known Issues Fixed

✅ Population no longer starts at 6 when using Default model
✅ Scenario is correctly set based on model (default/flood/creation)
✅ Each user gets their own copy of all models
✅ Parameters can be saved and loaded

## Next Steps for Frontend

1. **Model Dropdown**: Populate from `/api/models/list?username={user}`
2. **Load Parameters**: Call `/api/models/load` when model is selected
3. **Save Button**: POST to `/api/models/save` with modified parameters
4. **Reset Button**: POST to `/api/models/reset` to restore defaults
5. **Show Metadata**: Display base_model and last_modified from metadata

## Troubleshooting

**"User models not appearing"**
- Make sure you deleted the old `users/` directory
- Register a new user after rebuilding the executable
- Check `users/{username}/models/` directory was created

**"Population still starts at 6"**
- Check which model is selected in the dropdown
- Verify the model's parameters.csv has correct `start_pop_size`
- Make sure scenario is set to "default" not "flood"

**"Can't save/load parameters"**
- Verify the API endpoints are registered (check server.go)
- Check browser console for API errors
- Ensure metadata.json exists in model directory
