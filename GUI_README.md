# Drift GUI - User Guide

## Quick Start

### 1. Start the Web Server
```bash
drift.exe --web
```

### 2. Open Your Browser
Navigate to: `http://localhost:8080`

### 3. Login
- Select your username from the dropdown (default: alice)
- Leave password blank (optional)
- Click "Login"

### 4. Run a Simulation
1. Go to **Setup** tab
2. Configure parameters (or use defaults)
3. Click "Start Simulation"
4. Switch to **Progress** tab to watch it run
5. View results in **Graphs** tab when complete

## Features

### 🔐 Login Screen
- User selection from dropdown
- Optional password authentication
- Clean, modern interface

### ⚙️ Setup Tab
- **Basic Parameters**:
  - Model Name
  - End Year
  - Population sizes
  - Save interval
  - Number of runs

- **Tracking Options**:
  - Track DNA
  - Track on Map
  - Track Drift
  - Track Mutations

- **Actions**:
  - Start Simulation
  - Save Configuration (future)
  - Load Configuration (future)

### 📊 Progress Tab
- **Real-time Updates** (every 2 seconds):
  - Current status (queued/running/completed)
  - Current generation
  - Total generations
  - Population size

- **Progress Bar**:
  - Visual percentage completion
  - Smooth animations

- **Controls**:
  - Cancel simulation
  - Refresh progress manually

### 📈 Graphs Tab
- **Genetic Graph**: genetic.gif
- **Genealogical Graph**: genealogical.gif
- **Download Buttons**: Save results locally
- **Auto-refresh**: Images reload when simulation updates

### 🗺️ Maps Tab
- Displays map visualizations when "Track on Map" is enabled
- Shows info message when maps are disabled

## How It Works

### Login Flow
```
User selects username → API validates → Main app loads
```

### Simulation Flow
```
Setup parameters → Click Start → API queues job → drift.exe spawns
         ↓
Progress tab polls progress.json every 2 seconds
         ↓
Simulation completes → Graphs tab displays results
```

### File Structure
```
users/
└── alice/
    └── Default/
        ├── progress.json        # Real-time progress
        ├── Default_results.csv  # Results data
        ├── genetic.gif          # Genetic graph
        └── genealogical.gif     # Genealogical graph
```

## Technical Details

### API Endpoints Used
- `POST /api/login` - User authentication
- `GET /api/users` - Load user list
- `POST /api/simulation/start` - Start simulation
- `GET /api/progress/{user}/{model}` - Get progress
- `GET /api/results/{user}/{model}` - List result files
- `GET /api/results/{user}/{model}/{file}` - Download file
- `POST /api/simulation/cancel` - Cancel job

### Auto-refresh
- Progress tab polls every 2 seconds while simulation runs
- Stops polling when simulation completes
- Graphs reload with cache-busting timestamps

### State Management
JavaScript maintains global state:
- `currentUser` - Logged in username
- `currentModel` - Active model name
- `currentJobId` - Running job ID
- `progressInterval` - Polling interval

## User Experience

### Responsive Design
- Works on desktop, tablet, and mobile
- Adaptive grid layouts
- Touch-friendly buttons

### Visual Feedback
- Loading states
- Success/error messages
- Smooth transitions
- Progress animations

### Color Coding
- **Purple gradient**: Primary branding
- **Green**: Running status
- **Yellow**: Queued status
- **Blue**: Completed status
- **Red**: Error/cancel buttons

## Tips

### Best Practices
1. **Don't close browser** while simulation runs (progress will continue on server)
2. **Refresh Progress tab** to see latest updates
3. **Wait for completion** before downloading results
4. **Check Maps tab** only if track_map enabled

### Troubleshooting

**Can't login?**
- Check users/users.json exists
- Verify username exists in file
- Check browser console for errors

**Simulation won't start?**
- Ensure only one simulation per user at a time
- Check terminal for drift.exe errors
- Verify parameters are valid

**Progress not updating?**
- Click "Refresh" button
- Check if simulation is actually running (terminal)
- Verify progress.json exists in user directory

**Graphs not loading?**
- Wait for simulation to complete
- Check if .gif files exist in results directory
- Try refreshing the page

## Development

### File Locations
```
web/static/
├── index.html   # Main HTML structure
├── styles.css   # All styling
└── app.js       # Application logic
```

### Customization
- **Colors**: Edit CSS variables in styles.css
- **Poll interval**: Change `setInterval(checkProgress, 2000)` in app.js
- **Form fields**: Add inputs in index.html, update handleStartSimulation() in app.js

### Adding Parameters
1. Add input field in Setup tab HTML
2. Get value in `handleStartSimulation()`
3. Include in API request body
4. Backend already handles parameter passing to drift.exe

## Future Enhancements

### Planned Features
- [ ] Save/Load configuration presets
- [ ] Multiple model comparison
- [ ] Queue visibility (see other users' jobs)
- [ ] Admin dashboard
- [ ] Parameter validation
- [ ] Export all results as ZIP
- [ ] Dark mode toggle
- [ ] Email notifications on completion

### Nice-to-Have
- [ ] Real-time log streaming
- [ ] Interactive parameter tooltips
- [ ] Chart.js integration for additional graphs
- [ ] Session persistence (remember login)
- [ ] Drag-and-drop configuration upload

## Security Notes

### Current Implementation
- Basic username/password authentication
- No session tokens (stateless)
- No HTTPS (local only)
- Path traversal protection for file downloads

### For Production
- Add bcrypt password hashing
- Implement JWT tokens
- Enable HTTPS/TLS
- Add rate limiting
- Session management
- Admin role separation

## Browser Compatibility

✅ **Tested**:
- Chrome 90+
- Firefox 88+
- Edge 90+
- Safari 14+

⚠️ **Requires**:
- JavaScript enabled
- Modern browser (ES6+)
- Cookies enabled

## Performance

- **Lightweight**: ~50KB total (HTML + CSS + JS)
- **Fast**: Vanilla JavaScript, no frameworks
- **Efficient**: Polling stops when simulation completes
- **Responsive**: Instant UI updates

## Accessibility

- Semantic HTML
- Keyboard navigation
- Screen reader friendly
- High contrast text
- Focus indicators

---

**Built with ❤️ for the Drift genetic drift simulator**
