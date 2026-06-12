# Model Management System

## Overview

The simulation now uses a hierarchical model management system with:
- **Immutable base models** in `static/basemodels/`
- **User-specific model copies** in `users/{username}/models/{modelname}/`
- Full support for save/load/reset operations

## Directory Structure

```
static/basemodels/               # Immutable base templates
├── registry.csv                 # List of all base models
├── Default/
│   ├── parameters.csv
│   ├── actuarial_table.csv
│   └── chromosome_data.csv
├── QuickTest/
│   └── parameters.csv           # Inherits from Default
├── Flood/
│   └── parameters.csv           # Inherits from Default
├── Creation/
│   └── parameters.csv           # Inherits from Default
└── OoA/
    └── parameters.csv           # Inherits from Default

users/{username}/models/         # User's working copies
├── Default/
│   ├── metadata.json
│   ├── config/
│   │   └── parameters.csv
│   └── results/
├── QuickTest/
│   ├── metadata.json
│   ├── config/
│   └── results/
└── ...
```

## Base Models

### Default (Standard Human Population)
- Population: 1000
- Years: 100 (for testing, change to 4500 for full run)
- Scenario: default
- Uses age-distributed initialization

### QuickTest
- Population: 100
- Years: 500
- Scenario: default
- For rapid testing

### Flood (Biblical Flood)
- Population: 6 (3 couples)
- Years: 4500
- Scenario: flood
- Bottleneck at year 1656-1657

### Creation (Adam & Eve)
- Population: 2 (1 couple)
- Years: 4500
- Scenario: creation
- Extended lifespans (900 years)

### OoA (Out of Africa)
- Population: 10,000
- Years: 10,000
- Scenario: default
- For demographic modeling

## Parameter Inheritance

Models can inherit from each other using `_inherit` directive:

```csv
parameter,value
_inherit,Default
start_pop_size,100
end_year,500
```

This loads all parameters from Default, then overrides the specified values.

## Command Line Usage

```bash
# Load user model
drift -username=Rob -model=Default

# Load base model directly (temporary, not saved)
drift -base-model=Default

# Legacy: load from config directory
drift -config-root=static
```

## Web API

### List User Models
```
GET /api/models/list?username=Rob
```

### Load Model Parameters
```
GET /api/models/load?username=Rob&model_name=Default
```

### Save Model Parameters
```
POST /api/models/save
{
  "username": "Rob",
  "model_name": "Default",
  "parameters": {
    "end_year": 1000,
    "start_pop_size": 500
  }
}
```

### Reset Model to Base Defaults
```
POST /api/models/reset
{
  "username": "Rob",
  "model_name": "Default"
}
```

## User Initialization

When a new user is created, all base models are automatically copied to their `models/` directory. This ensures every user has access to all standard configurations.

## Scenarios

Different scenarios use different population initialization:

- **default**: Age-distributed population based on actuarial table
- **flood**: 3 pre-married couples with extended lifespans
- **creation**: 1 couple (Adam & Eve) with 900-year lifespans

## Key Files

- `pkg/models/manager.go` - Core model management logic
- `pkg/models/metadata.go` - Model metadata structures
- `pkg/models/user.go` - User initialization
- `pkg/config/modelmanager.go` - Config layer integration
- `pkg/webserver/handlers_models.go` - Web API handlers
- `static/basemodels/registry.csv` - Base model registry

## Benefits

✅ Separation of base templates from user data
✅ Easy reset to defaults
✅ Multiple models per user
✅ Version control friendly (CSV format)
✅ Parameter inheritance reduces duplication
✅ Clear provenance (metadata tracks base model)
✅ Web GUI support for save/load/reset

## Migration Notes

The old system stored models directly in `users/{username}/`. The new system uses `users/{username}/models/{modelname}/`. Both systems can coexist during transition.
