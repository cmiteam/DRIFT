# Plotting System Updates

## Changes Made

### 1. Fixed CSV Column Name Mapping
- **Problem**: The plotting code was looking for column names like `numinds`, `marriages`, etc., but the actual CSV uses `n`, `marrs`, `births`, `randDs`, `cullDs`, etc.
- **Solution**: Updated `METRIC_DEFINITIONS` to use the actual CSV column names from the simulation output

### 2. Removed Dependency on Setup Checkboxes
- **Problem**: Charts would only appear if you selected the corresponding checkboxes in Setup > Plot Options before running
- **Solution**:
  - Charts can now display ANY data from the CSV, regardless of setup selections
  - Setup checkboxes are now only used to suggest additional optional charts
  - All data is always available for custom charts

### 3. Always Create Standard Default Chart
- **New Feature**: Every simulation now gets a "Population Overview" chart with:
  - Population Size (n)
  - Births
  - Random Deaths (randDs)
  - Culled Deaths (cullDs)
- This chart is ALWAYS created, regardless of Plot Options selections

### 4. Optional Charts Based on Plot Options
- If you check any Ancestry metrics (Y/mt/genealogical/genetic descendants) in Setup, an "Ancestry Tracking" chart is created
- If you check any Genetic metrics (heterozygosity, mutations, fitness) in Setup, a "Genetic Metrics" chart is created
- These charts only include the specific metrics you selected

### 5. Enhanced Custom Chart Dialog
- **Add Chart** button now shows a proper dialog with:
  - Text input for chart name
  - Checkboxes for ALL available metrics from the CSV (not just the ones selected in setup)
  - Pretty formatting with colors matching the plot lines
  - Limit of 1-6 metrics per chart for readability
  - Shows metric display names (e.g., "Population Size" instead of "n")

### 6. Fixed X-Axis
- Changed X-axis from "Generation" to "Year" to match CSV data
- Uses the `year` column from the CSV

### 7. Expanded Metric Support
Added support for many more CSV columns:
- `nCents` - Centromeres
- `nAlleles` - Alleles
- `nBlocks` - Blocks
- `PercSeedGenoRet` - Seed Genome %
- `AvSeedGenoCov` - Seed Coverage
- `HomMin` - Homozygosity Min
- `HomMaj` - Homozygosity Maj
- And all the original metrics (n, births, deaths, ancestry, genetic metrics)

## User Workflow

### Before a Run
1. Go to Setup tab
2. Configure your simulation parameters
3. *Optionally* select metrics in Plot Options to get additional default charts
4. Start simulation

### During/After a Run
1. Go to Plots tab
2. See the standard "Population Overview" chart automatically
3. If you selected Plot Options, see additional charts for those metrics
4. Click **"+ Add Chart"** to create custom charts with ANY metrics
5. Toggle metrics on/off with checkboxes in chart headers
6. Click ⚙ to toggle Linear/Log scale
7. Click × to remove charts you don't need
8. Charts update every 3 seconds automatically

## Technical Details

### CSV Column Mapping
The system now correctly maps between:
- **Setup checkbox IDs** (e.g., `numinds`, `marriages`) used in the HTML form
- **CSV column names** (e.g., `n`, `marrs`) used in the actual data file
- **Display names** (e.g., "Population Size", "Marriages") shown to users

### Chart Structure
```javascript
{
    id: 'chart-id',
    title: 'Chart Title',
    metrics: ['n', 'births', 'randDs'],  // CSV column names
    yAxisType: 'linear' or 'log',
    isStandard: true/false  // Standard chart can't be deleted
}
```

### Data Flow
1. Simulation writes CSV with columns: `run`, `year`, `n`, `marrs`, `births`, etc.
2. Backend serves incremental data via `/api/plot-data/`
3. Frontend parses CSV and stores in `plotState.allData` keyed by column name
4. Charts reference metrics by CSV column name (e.g., `'n'`, `'births'`)
5. Display names are looked up from `METRIC_DEFINITIONS`

## Benefits

1. **No more empty plots** - Standard chart always shows population data
2. **Plot anything** - Not limited by what you selected in setup
3. **Flexible** - Add/remove charts during or after simulation
4. **User-friendly** - Clear dialogs and labels for custom charts
5. **Comprehensive** - Access to all CSV columns, not just a subset
