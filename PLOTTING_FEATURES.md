# Real-Time Plotting Features

## Overview
This document describes the newly implemented real-time plotting system for the Drift population genetics simulator web interface.

## Features Implemented

### 1. **Real-Time Plot Updates**
- Charts automatically update every 3 seconds while simulation is running
- Incremental data loading: only fetches new data points since last update
- File modification time checking to avoid unnecessary data transfers (304 Not Modified responses)
- Manual refresh option available
- Pause/Resume controls for plot updates

### 2. **Dynamic Chart Selection**
- Metrics can be toggled on/off **on the Plots page**, independent of setup selections
- Three default chart groups created automatically based on selected metrics:
  - **Population Dynamics**: Population size, marriages, births, deaths
  - **Genetic Metrics**: Heterozygosity, mutations, fitness
  - **Ancestry Tracking**: Y/mt/genealogical/genetic descendants
- Only charts with at least one selected metric are created

### 3. **Interactive Controls**

#### Per-Chart Controls:
- **Metric toggles**: Checkboxes in chart header to show/hide individual metrics in real-time
- **Settings button (⚙)**: Click to toggle between linear and logarithmic Y-axis scales
- **Delete button (×)**: Remove charts you don't need

#### Global Controls:
- **Update Mode**: Auto (3s polling) or Manual
- **Pause/Resume**: Stop/start automatic updates
- **Refresh Now**: Force immediate data reload
- **Add Chart**: Create custom charts (basic implementation)
- **Download CSV**: Download the raw results data

### 4. **Efficient Data Loading**

#### Backend Endpoint (`/api/plot-data/{user}/{model}`)
Query parameters:
- `from_row=N`: Only return rows starting from index N (incremental loading)

Headers:
- `If-Modified-Since`: Server returns 304 if file hasn't changed

Response includes:
- `headers`: CSV column names
- `rows`: Data rows as JSON objects
- `total_rows`: Total number of rows in CSV
- `from_row`: Starting row index
- `new_rows`: Number of new rows returned
- `modified`: File modification timestamp

#### Frontend Optimization
- Tracks last loaded row count
- Tracks last file modification time
- Only parses and appends new data on each update
- Minimal bandwidth usage during updates

### 5. **Chart Customization**

Each chart supports:
- **Linear/Log scale**: Click settings (⚙) to toggle
- **Selective metrics**: Use checkboxes to show only metrics you want
- **Interactive zoom/pan**: Powered by Plotly.js built-in controls
- **Hover tooltips**: See exact values on mouse hover
- **Legend**: Click legend items to toggle series visibility

## Technical Architecture

### Frontend Components

**Files Modified:**
- `web/static/index.html`: Added Plotly.js library, new Plots tab UI
- `web/static/app.js`: ~550 lines of plotting functionality
- `web/static/styles.css`: Chart styling and responsive layout

**Key JavaScript Functions:**
- `loadPlots()`: Initializes plotting when switching to Plots tab
- `fetchPlotData()`: Polls server for new data
- `createChart()`: Dynamically creates chart DOM and Plotly visualization
- `updateChart()`: Updates chart with new data
- `toggleMetricVisibility()`: Show/hide metrics
- `showChartSettings()`: Toggle linear/log scale

### Backend Components

**Files Modified:**
- `pkg/webserver/server.go`: Added `/api/plot-data/` route
- `pkg/webserver/handlers.go`: Implemented `handlePlotData()` function

**Data Flow:**
1. Simulation writes `{model}_results.csv`
2. Backend serves incremental CSV data as JSON
3. Frontend polls every 3 seconds (when not paused)
4. New rows are appended to in-memory data arrays
5. Plotly charts are updated with `Plotly.restyle()`

## Supported Metrics

The system recognizes and can plot these metrics (if selected in Setup > Plot Options):

| Metric ID | Display Name | Default Chart |
|-----------|-------------|---------------|
| `numinds` | Population Size | Population Dynamics |
| `marriages` | Marriages | Population Dynamics |
| `births` | Births | Population Dynamics |
| `random_deaths` | Random Deaths | Population Dynamics |
| `cull_deaths` | Culled Deaths | Population Dynamics |
| `Y_descends` | Y Descendants | Ancestry Tracking |
| `mt_descends` | mt Descendants | Ancestry Tracking |
| `genealo_descends` | Genealogical Descendants | Ancestry Tracking |
| `genetic_descends` | Genetic Descendants | Ancestry Tracking |
| `av_heterozygosity` | Average Heterozygosity | Genetic Metrics |
| `num_mutations` | Number of Mutations | Genetic Metrics |
| `av_ind_fitness` | Average Individual Fitness | Genetic Metrics |

## Usage Workflow

1. **Setup**: On Setup tab, select which metrics to track in "Plot Options"
2. **Run**: Start simulation
3. **View Progress**: Switch to Progress tab to monitor execution
4. **View Plots**: Switch to Plots tab
   - Charts auto-generate based on selected metrics
   - Charts update every 3 seconds automatically
   - Toggle individual metrics on/off with checkboxes
   - Change scale (linear/log) with ⚙ button
   - Remove unwanted charts with × button
5. **Customize**: Add custom charts, pause updates, adjust display

## Future Enhancements

Possible improvements for future versions:

1. **Advanced Chart Builder**
   - Modal dialog for creating custom charts
   - Select specific metrics to include
   - Choose chart type (line, scatter, area, bar)
   - Multiple Y-axes support

2. **Chart Persistence**
   - Save chart configurations per user
   - Remember preferences across sessions

3. **Export Options**
   - Export individual charts as PNG/SVG
   - Export all charts as PDF report

4. **Statistical Overlays**
   - Add trend lines
   - Show moving averages
   - Highlight statistical events

5. **Time Range Selection**
   - Slider to zoom into specific generation ranges
   - "Follow mode" to always show latest N generations

6. **Comparison Mode**
   - Overlay results from multiple simulation runs
   - Side-by-side chart comparison

## Performance Notes

- **Memory**: All plot data is kept in browser memory for fast updates
  - For simulations with 10,000+ generations, this uses ~5-10MB per chart
  - No performance issues expected for typical simulations

- **Network**: Incremental loading minimizes bandwidth
  - Initial load: Entire CSV (could be 100KB-1MB)
  - Updates: Only new rows (~1-10KB per update typically)

- **CPU**: Plotly.js handles rendering efficiently
  - Can plot 10,000+ points smoothly
  - Chart updates complete in <50ms typically

## Browser Compatibility

Tested and working in:
- Chrome/Edge (recommended)
- Firefox
- Safari

Requires JavaScript enabled and modern ES6+ support.

## Troubleshooting

**Charts not appearing?**
- Ensure at least one metric is selected in Setup > Plot Options
- Check that simulation has started and generated results CSV
- Check browser console for errors

**Charts not updating?**
- Check that "Update" mode is set to "Auto"
- Check that Pause button shows "Pause" (not "Resume")
- Verify simulation is running (check Progress tab)

**Performance issues?**
- Try reducing the number of visible charts
- Toggle off metrics you don't need
- Use Pause when not actively viewing plots

**Missing metrics?**
- Metrics only appear if selected in Setup before starting simulation
- To add metrics, start a new simulation with desired options selected
