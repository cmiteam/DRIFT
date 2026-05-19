// Global State
const state = {
    currentUser: null,
    currentModel: 'Default',
    currentJobId: null,
    progressInterval: null,
    loadedParameters: {}  // Track parameters as they were loaded
};

// API Base URL
const API_BASE = '/api';

// Parameter metadata definitions - specifies how parameters should be rendered
// Parameters not in this list will use default rendering based on their value type
const PARAMETER_METADATA = {
    // Basic parameters (shown in Basic tab)
    'num_runs': { tab: 'basic', label: 'Number of Runs', type: 'int', min: 1 },
    'start_pop_size': { tab: 'basic', label: 'Start Population Size', type: 'int', min: 1 },
    'max_pop_size': { tab: 'basic', label: 'Max Population Size', type: 'int', min: 1 },
    'end_year': { tab: 'basic', label: 'End Year', type: 'int', min: 1 },
    'save_interval': { tab: 'basic', label: 'Save Interval', type: 'int', min: 1 },
    'model_name': { tab: 'basic', label: 'Model Name', type: 'string' },
    'base_model': { tab: 'basic', label: 'Base Model', type: 'dropdown' },

    // DNA parameters (shown in DNA tab)
    'track_DNA': { tab: 'dna', label: 'Track DNA', type: 'checkbox' },
    'seed_year': { tab: 'dna', label: 'Seed Year', type: 'int', min: 0 },
    'seed_style': { tab: 'dna', label: 'Seed Style', type: 'dropdown', options: [
        { value: 0, label: 'Single' },
        { value: 1, label: 'Population' },
        { value: 2, label: 'Max Heterozygosity' }
    ]},
    'init_heterozygosity': { tab: 'dna', label: 'Initial Heterozygosity', type: 'float', min: 0, max: 1, step: 0.01 },
    'multiplier': { tab: 'dna', label: 'Multiplier', type: 'int', min: 1 },
    'genome_map': { tab: 'dna', label: 'Genome Map', type: 'checkbox' },
    'every_genome_map': { tab: 'dna', label: 'All Genome Maps', type: 'checkbox' },

    // Mutation parameters (shown in Mutation tab)
    'track_mutations': { tab: 'mutation', label: 'Track Mutations', type: 'checkbox' },
    'mu': { tab: 'mutation', label: 'Mutation Rate', type: 'float', min: 0, step: 0.1 },
    'selection': { tab: 'mutation', label: 'Selection', type: 'dropdown', options: [
        { value: 0, label: 'None' },
        { value: 1, label: 'Active' }
    ]},
    'f_neutral': { tab: 'mutation', label: 'f(Neutral)', type: 'float', min: 0, max: 1, step: 0.0001 },
    'f_beneficial': { tab: 'mutation', label: 'f(Beneficial)', type: 'float', min: 0, max: 1, step: 0.0001 },
    'mu_scale_factor': { tab: 'mutation', label: 'Mutation Scale Factor', type: 'int', min: 1 },
    'shape': { tab: 'mutation', label: 'Weibull Shape', type: 'float', min: 0, step: 0.01 },
    'scale': { tab: 'mutation', label: 'Weibull Scale', type: 'float', min: 0, step: 0.01 },
    'Weibull_adj': { tab: 'mutation', label: 'Weibull Adjustment', type: 'int', min: 1 },
    'mutation_map': { tab: 'mutation', label: 'Mutation Map', type: 'checkbox' },
    'mutation_hist': { tab: 'mutation', label: 'Mutation Histogram', type: 'checkbox' },

    // Map parameters (shown in Maps tab)
    'track_map': { tab: 'maps', label: 'Track on Map', type: 'checkbox' },
    'map_name': { tab: 'maps', label: 'Map Name', type: 'dropdown' },
    'scaling': { tab: 'maps', label: 'Scaling', type: 'float', min: 0.1, step: 0.1 },

    // Tracking parameters (shown in Tracking tab)
    'track_drift': { tab: 'tracking', label: 'Track Drift', type: 'checkbox' },
    'track_dead': { tab: 'tracking', label: 'Track Dead', type: 'checkbox' },
    'track_SFS': { tab: 'tracking', label: 'Site Frequency Spectrum', type: 'checkbox' },
    'save_detailed_SFS': { tab: 'tracking', label: 'Save Detailed SFS', type: 'checkbox' },
    'track_LD': { tab: 'tracking', label: 'Linkage Disequilibrium', type: 'checkbox' },
    'track_coalescence': { tab: 'tracking', label: 'Track Coalescence', type: 'checkbox' },

    // Other parameters (shown in Other Parameters tab)
    'max_growth_rate': { tab: 'other', label: 'Max Growth Rate', type: 'float', min: 0, step: 0.00001, group: 'Population Dynamics' },
    'bottleneck_start': { tab: 'other', label: 'Bottleneck Start Year', type: 'int', description: '-1 to disable', group: 'Bottleneck' },
    'bottleneck_end': { tab: 'other', label: 'Bottleneck End Year', type: 'int', min: 0, group: 'Bottleneck' },
    'bottleneck_size': { tab: 'other', label: 'Bottleneck Size', type: 'int', min: 1, group: 'Bottleneck' },
    'lifespan': { tab: 'other', label: 'Lifespan', type: 'int', min: 1, group: 'Life Stage' },
    'min_lifespan': { tab: 'other', label: 'Min Lifespan', type: 'int', min: 1, group: 'Life Stage' },
    'lifespan_drop': { tab: 'other', label: 'Lifespan Drop', type: 'float', min: 0, max: 1, step: 0.01, group: 'Life Stage' },
    'maturity': { tab: 'other', label: 'Maturity Age', type: 'int', min: 0, group: 'Life Stage' },
    'menopause': { tab: 'other', label: 'Menopause Ratio', type: 'float', min: 0, max: 1, step: 0.01, group: 'Reproduction' },
    'birth_prob': { tab: 'other', label: 'Birth Probability', type: 'float', min: 0, step: 0.1, group: 'Reproduction' },
    'spacing': { tab: 'other', label: 'Birth Spacing', type: 'int', min: 1, group: 'Reproduction' },
    'mating_style': { tab: 'other', label: 'Mating Style', type: 'dropdown', options: [
        { value: 0, label: 'Random' },
        { value: 1, label: 'Distance-based' },
        { value: 2, label: 'Age and Distance' }
    ], group: 'Mating & Geography' },
    'max_mating_distance': { tab: 'other', label: 'Max Mating Distance', type: 'int', min: 0, group: 'Mating & Geography' },
    'residence_style': { tab: 'other', label: 'Residence Style', type: 'int', min: 0, group: 'Mating & Geography' },
    'wander': { tab: 'other', label: 'Wander Distance', type: 'int', min: 0, group: 'Mating & Geography' },
    'max_breeding_inds': { tab: 'other', label: 'Max Breeding Individuals', type: 'int', description: '-1 for unlimited', group: 'Mating & Geography' },
    'sex_chrom_index': { tab: 'other', label: 'Sex Chromosome Index', type: 'int', min: 1, group: 'Genetics' },
    'animation_save_interval': { tab: 'other', label: 'Animation Save Interval', type: 'int', min: 1, group: 'Save Settings' }
};

// Parameters that should never be shown in Other Parameters (already handled elsewhere)
const EXCLUDED_FROM_OTHER = new Set([
    'username', 'scenario', 'base_model', 'model_name',
    // Basic tab
    'num_runs', 'start_pop_size', 'max_pop_size', 'end_year', 'save_interval',
    // DNA tab
    'track_DNA', 'seed_year', 'seed_style', 'init_heterozygosity', 'multiplier', 'genome_map', 'every_genome_map',
    // Mutation tab
    'track_mutations', 'mu', 'selection', 'f_neutral', 'f_beneficial', 'mu_scale_factor', 'shape', 'scale', 'Weibull_adj', 'mutation_map', 'mutation_hist',
    // Maps tab
    'track_map', 'map_name', 'scaling',
    // Tracking tab
    'track_drift', 'track_dead', 'track_SFS', 'save_detailed_SFS', 'track_LD', 'track_coalescence'
]);

// Initialize App
document.addEventListener('DOMContentLoaded', () => {
    initializeApp();
});

async function initializeApp() {
    // Setup event listeners
    setupEventListeners();
    setupFormChangeListeners();
}

// Toggle between login and registration forms
function showRegistrationForm() {
    document.getElementById('loginForm').style.display = 'none';
    document.getElementById('showRegisterBtn').style.display = 'none';
    document.querySelector('.divider').style.display = 'none';
    document.getElementById('registerForm').style.display = 'block';
    document.getElementById('loginError').classList.remove('show');
    document.getElementById('registerSuccess').classList.remove('show');
}

function showLoginForm() {
    document.getElementById('loginForm').style.display = 'block';
    document.getElementById('showRegisterBtn').style.display = 'block';
    document.querySelector('.divider').style.display = 'flex';
    document.getElementById('registerForm').style.display = 'none';
    document.getElementById('loginError').classList.remove('show');
    document.getElementById('registerSuccess').classList.remove('show');
    document.getElementById('registerForm').reset();
}

// Setup Event Listeners
function setupEventListeners() {
    // Login form
    document.getElementById('loginForm').addEventListener('submit', handleLogin);

    // Registration form toggle
    document.getElementById('showRegisterBtn').addEventListener('click', showRegistrationForm);
    document.getElementById('cancelRegisterBtn').addEventListener('click', showLoginForm);
    document.getElementById('registerForm').addEventListener('submit', handleRegistration);

    // Logout button
    document.getElementById('logoutBtn').addEventListener('click', handleLogout);

    // Tab buttons
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.addEventListener('click', () => switchTab(btn.dataset.tab));
    });

    // Parameter tab buttons
    document.querySelectorAll('.param-tab-btn').forEach(btn => {
        btn.addEventListener('click', () => switchParamTab(btn.dataset.paramTab));
    });

    // Enable/disable dependent content
    setupDependentContent('track_map', 'mapsParamsContent');
    setupDependentContent('track_DNA', 'dnaParamsContent');
    setupDependentContent('track_mutations', 'mutationParamsContent');

    // Setup form
    document.getElementById('setupForm').addEventListener('submit', handleStartSimulation);

    // Base model dropdown change handler
    document.getElementById('base_model')?.addEventListener('change', handleBaseModelChange);

    // Model management buttons
    document.getElementById('saveConfig')?.addEventListener('click', handleSaveModel);
    document.getElementById('resetToDefaults')?.addEventListener('click', handleResetToDefaults);
    document.getElementById('saveAsBaseModel')?.addEventListener('click', handleSaveAsNewModel);

    // Progress controls
    document.getElementById('cancelSim')?.addEventListener('click', handleCancelSimulation);
    document.getElementById('refreshProgress')?.addEventListener('click', checkProgress);

    // Download buttons
    document.querySelectorAll('.download-btn').forEach(btn => {
        btn.addEventListener('click', () => downloadFile(btn.dataset.file));
    });

    document.getElementById('downloadCSV')?.addEventListener('click', () => {
        const modelName = state.currentModel || 'Default';
        downloadFile(`${modelName}_results.csv`);
    });
}

// Handle Login
async function handleLogin(e) {
    e.preventDefault();

    const username = document.getElementById('loginUsername').value;
    const password = document.getElementById('loginPassword').value;
    const errorDiv = document.getElementById('loginError');

    console.log('Login attempt:', username);

    try {
        const response = await fetch(`${API_BASE}/login`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ username, password })
        });

        console.log('Login response status:', response.status);
        const data = await response.json();
        console.log('Login response data:', data);

        if (data.success) {
            state.currentUser = username;
            showMainApp();
            checkForExistingSimulation();
        } else {
            errorDiv.textContent = data.message || 'Login failed';
            errorDiv.classList.add('show');
        }
    } catch (error) {
        let errorMessage = error.message;
        if (error.message === 'Failed to fetch') {
            errorMessage = 'Cannot connect to server. Please ensure the server is running (go run drift.go -web)';
        }
        errorDiv.textContent = 'Login error: ' + errorMessage;
        errorDiv.classList.add('show');
    }
}

// Handle Registration
async function handleRegistration(e) {
    e.preventDefault();

    const username = document.getElementById('registerUsername').value;
    const password = document.getElementById('registerPassword').value;
    const errorDiv = document.getElementById('loginError');
    const successDiv = document.getElementById('registerSuccess');

    // Clear previous messages
    errorDiv.classList.remove('show');
    successDiv.classList.remove('show');

    console.log('Registering user:', username);

    try {
        const response = await fetch(`${API_BASE}/register`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ username, password })
        });

        console.log('Response status:', response.status);
        const data = await response.json();
        console.log('Response data:', data);

        if (data.success) {
            // Check if approval is required
            if (data.requires_approval) {
                successDiv.textContent = 'Registration successful! Your account is pending administrator approval.';
                successDiv.classList.add('show');
            } else {
                successDiv.textContent = 'Registration successful! You can now log in.';
                successDiv.classList.add('show');
            }

            // Reset form and switch back to login after 2 seconds
            setTimeout(() => {
                showLoginForm();
                // Pre-fill username in login form
                document.getElementById('loginUsername').value = username;
            }, 2000);
        } else {
            errorDiv.textContent = data.message || 'Registration failed';
            errorDiv.classList.add('show');
        }
    } catch (error) {
        let errorMessage = error.message;
        if (error.message === 'Failed to fetch') {
            errorMessage = 'Cannot connect to server. Please ensure the server is running (go run drift.go -web)';
        }
        errorDiv.textContent = 'Registration error: ' + errorMessage;
        errorDiv.classList.add('show');
    }
}

// Handle Logout
function handleLogout() {
    state.currentUser = null;
    state.currentJobId = null;
    clearInterval(state.progressInterval);

    document.getElementById('loginScreen').classList.add('active');
    document.getElementById('mainApp').classList.remove('active');

    // Reset forms and show login screen
    document.getElementById('loginForm').reset();
    showLoginForm();
}

// Show Main App
function showMainApp() {
    console.log('showMainApp called, user:', state.currentUser);
    const loginScreen = document.getElementById('loginScreen');
    const mainApp = document.getElementById('mainApp');
    const currentUser = document.getElementById('currentUser');

    console.log('loginScreen:', loginScreen);
    console.log('mainApp:', mainApp);
    console.log('currentUser:', currentUser);

    loginScreen.classList.remove('active');
    mainApp.classList.add('active');
    currentUser.textContent = `Logged in as: ${state.currentUser}`;

    console.log('UI updated');

    // Load available models and maps
    loadBaseModels();
    loadAvailableModels();
    loadAvailableMaps();
}

// Load Base Models (templates for new models)
async function loadBaseModels() {
    try {
        const response = await fetch(`${API_BASE}/basemodels/list`);
        if (response.ok) {
            const data = await response.json();
            const baseModelSelect = document.getElementById('base_model');

            if (data.models && data.models.length > 0) {
                // Clear existing options
                baseModelSelect.innerHTML = '';

                // Add base models to dropdown
                data.models.forEach(model => {
                    const option = document.createElement('option');
                    option.value = model.id;
                    option.textContent = model.name;
                    option.title = model.description;
                    baseModelSelect.appendChild(option);
                });

                // Auto-load the first model - try user's saved version first
                const firstModelId = data.models[0].id;
                let loaded = false;

                if (state.currentUser) {
                    try {
                        const userResponse = await fetch(`${API_BASE}/models/load?username=${state.currentUser}&model_name=${firstModelId}`);
                        if (userResponse.ok) {
                            const userData = await userResponse.json();
                            if (userData.success && userData.parameters) {
                                populateFormFromParameters(userData.parameters);
                                document.getElementById('model_name').value = firstModelId;
                                state.currentModel = firstModelId;
                                loaded = true;
                            }
                        }
                    } catch (error) {
                        console.log(`No saved user model for ${firstModelId}`);
                    }
                }

                // Fall back to base model defaults
                if (!loaded) {
                    const paramsResponse = await fetch(`${API_BASE}/basemodels/parameters?id=${firstModelId}`);
                    if (paramsResponse.ok) {
                        const paramsData = await paramsResponse.json();
                        if (paramsData.success && paramsData.parameters) {
                            populateFormFromParameters(paramsData.parameters);
                            document.getElementById('model_name').value = firstModelId;
                            state.currentModel = firstModelId;
                        }
                    }
                }
            }
        }
    } catch (error) {
        console.error('Failed to load base models:', error);
    }
}

// Handle Base Model Selection Change
async function handleBaseModelChange(e) {
    const baseModelId = e.target.value;
    if (!baseModelId) return;

    // Try to load user's saved model first, fall back to base model defaults
    if (state.currentUser) {
        try {
            const userResponse = await fetch(`${API_BASE}/models/load?username=${state.currentUser}&model_name=${baseModelId}`);
            if (userResponse.ok) {
                const userData = await userResponse.json();
                if (userData.success && userData.parameters) {
                    populateFormFromParameters(userData.parameters);
                    document.getElementById('model_name').value = baseModelId;
                    state.currentModel = baseModelId;
                    console.log(`Loaded user's saved model: ${baseModelId}`);
                    return;
                }
            }
        } catch (error) {
            console.log(`No saved user model for ${baseModelId}, loading base defaults`);
        }
    }

    // Fall back to base model defaults
    try {
        const response = await fetch(`${API_BASE}/basemodels/parameters?id=${baseModelId}`);
        if (response.ok) {
            const data = await response.json();
            if (data.success && data.parameters) {
                populateFormFromParameters(data.parameters);
                document.getElementById('model_name').value = baseModelId;
                state.currentModel = baseModelId;
            }
        } else {
            console.error('Failed to load base model parameters');
        }
    } catch (error) {
        console.error('Failed to load base model parameters:', error);
    }
}

// Load Available Models
async function loadAvailableModels() {
    if (!state.currentUser) {
        console.error('No user logged in');
        return;
    }

    try {
        const response = await fetch(`${API_BASE}/models/list?username=${state.currentUser}`);
        if (response.ok) {
            const data = await response.json();

            if (data.success && data.models && data.models.length > 0) {
                // Load the first model's parameters by default
                if (data.models.length > 0) {
                    state.currentModel = data.models[0].id;
                    await loadModelParameters(data.models[0].id);
                }
            }
        }
    } catch (error) {
        console.error('Failed to load available models:', error);
    }
}

// Load parameters for a specific model
async function loadModelParameters(modelName) {
    if (!state.currentUser) {
        console.error('No user logged in');
        return;
    }

    try {
        const response = await fetch(`${API_BASE}/models/load?username=${state.currentUser}&model_name=${modelName}`);
        if (response.ok) {
            const data = await response.json();
            if (data.success && data.parameters) {
                // Populate form fields with loaded parameters
                populateFormFromParameters(data.parameters);
                console.log(`Loaded parameters for model: ${modelName}`);
            }
        }
    } catch (error) {
        console.error('Failed to load model parameters:', error);
    }
}

// Populate form fields from parameter object
function populateFormFromParameters(params) {
    // Store loaded parameters for change detection
    state.loadedParameters = { ...params };

    // First, render dynamic "Other Parameters" for any params not in existing tabs
    renderOtherParameters(params);

    // Iterate through all parameters and set form values
    for (const [key, value] of Object.entries(params)) {
        const element = document.getElementById(key);
        if (element) {
            if (element.type === 'checkbox') {
                element.checked = value == 1 || value === true;
            } else {
                element.value = value;
            }
        }
    }

    // Trigger dependent content updates
    updateDependentContent('track_map', 'mapsParamsContent');
    updateDependentContent('track_DNA', 'dnaParamsContent');
    updateDependentContent('track_mutations', 'mutationParamsContent');

    // Reset button states (no changes yet)
    updateChangeButtons();
}

// Render dynamic parameters in the "Other Parameters" tab
function renderOtherParameters(params) {
    const container = document.getElementById('otherParamsContainer');
    if (!container) return;

    // Clear existing content
    container.innerHTML = '';

    // Group parameters by their group property
    const paramsByGroup = {};
    const ungroupedParams = [];

    for (const [key, value] of Object.entries(params)) {
        // Skip parameters that are already handled in other tabs
        if (EXCLUDED_FROM_OTHER.has(key)) continue;

        const metadata = PARAMETER_METADATA[key];

        // If metadata says it belongs to a different tab, skip it
        if (metadata && metadata.tab && metadata.tab !== 'other') continue;

        // Group the parameter
        const group = metadata?.group || 'Other';
        if (group === 'Other') {
            ungroupedParams.push({ key, value, metadata });
        } else {
            if (!paramsByGroup[group]) {
                paramsByGroup[group] = [];
            }
            paramsByGroup[group].push({ key, value, metadata });
        }
    }

    // Render grouped parameters
    const groupOrder = ['Population Dynamics', 'Bottleneck', 'Life Stage', 'Reproduction', 'Mating & Geography', 'Genetics', 'Save Settings'];

    for (const groupName of groupOrder) {
        const groupParams = paramsByGroup[groupName];
        if (!groupParams || groupParams.length === 0) continue;

        const groupDiv = document.createElement('div');
        groupDiv.className = 'param-group';

        const groupHeader = document.createElement('h4');
        groupHeader.textContent = groupName;
        groupDiv.appendChild(groupHeader);

        const formRow = document.createElement('div');
        formRow.className = 'form-row dynamic-params';

        for (const param of groupParams) {
            const fieldDiv = createParameterField(param.key, param.value, param.metadata);
            formRow.appendChild(fieldDiv);
        }

        groupDiv.appendChild(formRow);
        container.appendChild(groupDiv);
    }

    // Render any remaining ungrouped parameters (parameters not in PARAMETER_METADATA)
    if (ungroupedParams.length > 0) {
        const ungroupedDiv = document.createElement('div');
        ungroupedDiv.className = 'param-group';

        const ungroupedHeader = document.createElement('h4');
        ungroupedHeader.textContent = 'Additional Settings';
        ungroupedDiv.appendChild(ungroupedHeader);

        const formRow = document.createElement('div');
        formRow.className = 'form-row dynamic-params';

        for (const param of ungroupedParams) {
            const fieldDiv = createParameterField(param.key, param.value, param.metadata);
            formRow.appendChild(fieldDiv);
        }

        ungroupedDiv.appendChild(formRow);
        container.appendChild(ungroupedDiv);
    }

    // Add change listeners to newly created elements
    container.querySelectorAll('input, select').forEach(element => {
        element.addEventListener('change', updateChangeButtons);
        element.addEventListener('input', updateChangeButtons);
    });
}

// Create a form field for a parameter based on its type and metadata
function createParameterField(key, value, metadata) {
    const fieldDiv = document.createElement('div');
    fieldDiv.className = 'form-group';

    const label = document.createElement('label');
    label.setAttribute('for', key);
    label.textContent = metadata?.label || formatParameterName(key);

    // Add description as tooltip if provided
    if (metadata?.description) {
        label.title = metadata.description;
        label.textContent += ' ⓘ';
    }

    fieldDiv.appendChild(label);

    // Determine field type based on metadata or value
    const fieldType = determineFieldType(key, value, metadata);
    let input;

    switch (fieldType) {
        case 'checkbox':
            input = document.createElement('input');
            input.type = 'checkbox';
            input.id = key;
            input.value = '1';
            input.checked = value == 1 || value === true;
            // Wrap checkbox with its label
            const checkboxLabel = document.createElement('label');
            checkboxLabel.className = 'checkbox-wrapper';
            checkboxLabel.appendChild(input);
            checkboxLabel.appendChild(document.createTextNode(' ' + (metadata?.label || formatParameterName(key))));
            fieldDiv.innerHTML = '';
            fieldDiv.appendChild(checkboxLabel);
            break;

        case 'dropdown':
            input = document.createElement('select');
            input.id = key;
            if (metadata?.options) {
                for (const opt of metadata.options) {
                    const option = document.createElement('option');
                    option.value = opt.value;
                    option.textContent = opt.label;
                    if (opt.value == value) {
                        option.selected = true;
                    }
                    input.appendChild(option);
                }
            }
            fieldDiv.appendChild(input);
            break;

        case 'float':
            input = document.createElement('input');
            input.type = 'number';
            input.id = key;
            input.value = value;
            input.step = metadata?.step || 'any';
            if (metadata?.min !== undefined) input.min = metadata.min;
            if (metadata?.max !== undefined) input.max = metadata.max;
            fieldDiv.appendChild(input);
            break;

        case 'int':
        default:
            input = document.createElement('input');
            input.type = 'number';
            input.id = key;
            input.value = value;
            input.step = '1';
            if (metadata?.min !== undefined) input.min = metadata.min;
            if (metadata?.max !== undefined) input.max = metadata.max;
            fieldDiv.appendChild(input);
            break;
    }

    return fieldDiv;
}

// Determine field type based on metadata and value
function determineFieldType(key, value, metadata) {
    // If metadata specifies type, use it
    if (metadata?.type) {
        return metadata.type;
    }

    // Infer from key name
    if (key.startsWith('track_') || key.endsWith('_map') || key.endsWith('_hist')) {
        return 'checkbox';
    }

    // Infer from value
    if (typeof value === 'boolean' || value === 0 || value === 1) {
        // Check if it's likely a boolean based on key name patterns
        if (key.includes('track') || key.includes('save_detailed') || key.includes('enable')) {
            return 'checkbox';
        }
    }

    // Check if it's a float
    if (typeof value === 'number' && !Number.isInteger(value)) {
        return 'float';
    }

    return 'int';
}

// Format parameter name for display (convert snake_case to Title Case)
function formatParameterName(key) {
    return key
        .split('_')
        .map(word => word.charAt(0).toUpperCase() + word.slice(1))
        .join(' ');
}

// Check if form has changes compared to loaded parameters
function hasFormChanges() {
    for (const [key, loadedValue] of Object.entries(state.loadedParameters)) {
        const element = document.getElementById(key);
        if (element) {
            let currentValue;
            if (element.type === 'checkbox') {
                currentValue = element.checked ? 1 : 0;
            } else {
                currentValue = element.value;
            }
            // Compare as strings to handle type differences
            if (String(currentValue) !== String(loadedValue)) {
                return true;
            }
        }
    }
    // Also check model_name field
    const modelNameEl = document.getElementById('model_name');
    const baseModelEl = document.getElementById('base_model');
    if (modelNameEl && baseModelEl && modelNameEl.value !== state.currentModel) {
        return true;
    }
    return false;
}

// Update Save and Reset button styles based on changes
function updateChangeButtons() {
    const saveBtn = document.getElementById('saveConfig');
    const resetBtn = document.getElementById('resetToDefaults');

    const hasChanges = hasFormChanges();

    if (saveBtn) {
        saveBtn.classList.toggle('btn-changed', hasChanges);
    }
    if (resetBtn) {
        resetBtn.classList.toggle('btn-changed', hasChanges);
    }
}

// Setup form change listeners
function setupFormChangeListeners() {
    const form = document.getElementById('setupForm');
    if (!form) return;

    // Listen for changes on all inputs, selects, and checkboxes
    form.querySelectorAll('input, select').forEach(element => {
        element.addEventListener('change', updateChangeButtons);
        element.addEventListener('input', updateChangeButtons);
    });
}

// Load Available Maps
async function loadAvailableMaps() {
    try {
        const response = await fetch(`${API_BASE}/maps/list`);
        if (response.ok) {
            const data = await response.json();
            const mapSelect = document.getElementById('map_name');

            if (data.maps && data.maps.length > 0) {
                // Clear existing options
                mapSelect.innerHTML = '';

                // Add maps to dropdown
                data.maps.forEach(map => {
                    const option = document.createElement('option');
                    option.value = map;
                    option.textContent = map.charAt(0).toUpperCase() + map.slice(1);
                    mapSelect.appendChild(option);
                });
            }
        }
    } catch (error) {
        console.error('Failed to load available maps:', error);
        // Keep default "sandbox" option if loading fails
    }
}

// Switch Tabs
function switchTab(tabName) {
    // Update tab buttons
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.classList.remove('active');
        if (btn.dataset.tab === tabName) {
            btn.classList.add('active');
        }
    });

    // Update tab panes
    document.querySelectorAll('.tab-pane').forEach(pane => {
        pane.classList.remove('active');
    });
    document.getElementById(`${tabName}Tab`).classList.add('active');

    // Clear status messages when switching away from setup
    if (tabName !== 'setup') {
        const statusDiv = document.getElementById('setupStatus');
        if (statusDiv) {
            statusDiv.textContent = '';
            statusDiv.className = 'status-message';
        }
    }

    // Refresh data when switching to certain tabs
    if (tabName === 'progress') {
        checkProgress();
    } else if (tabName === 'plots') {
        loadPlots();
    } else if (tabName === 'maps') {
        loadMaps();
    }
}

// Switch Parameter Tabs
function switchParamTab(tabName) {
    // Update parameter tab buttons
    document.querySelectorAll('.param-tab-btn').forEach(btn => {
        btn.classList.remove('active');
        if (btn.dataset.paramTab === tabName) {
            btn.classList.add('active');
        }
    });

    // Update parameter tab content
    document.querySelectorAll('.param-tab-content').forEach(content => {
        content.classList.remove('active');
    });

    // Map tab names to their corresponding content divs
    const tabMapping = {
        'basic': 'basicParams',
        'maps': 'mapsParams',
        'dna': 'dnaParams',
        'mutation': 'mutationParams',
        'tracking': 'trackingParams',
        'other': 'otherParams'
    };

    const contentId = tabMapping[tabName];
    if (contentId) {
        document.getElementById(contentId).classList.add('active');
    }
}

// Setup Dependent Content Enable/Disable
function setupDependentContent(checkboxId, contentId) {
    const checkbox = document.getElementById(checkboxId);
    const content = document.getElementById(contentId);

    if (!checkbox || !content) return;

    // Function to update content state
    const updateContentState = () => {
        if (checkbox.checked) {
            content.classList.remove('disabled');
            // Enable all inputs within the content
            content.querySelectorAll('input, select').forEach(input => {
                input.disabled = false;
            });
        } else {
            content.classList.add('disabled');
            // Disable all inputs within the content
            content.querySelectorAll('input, select').forEach(input => {
                input.disabled = true;
            });
        }
    };

    // Set initial state
    updateContentState();

    // Add event listener
    checkbox.addEventListener('change', updateContentState);
}

// Update dependent content state without setting up event listeners
function updateDependentContent(checkboxId, contentId) {
    const checkbox = document.getElementById(checkboxId);
    const content = document.getElementById(contentId);

    if (!checkbox || !content) return;

    if (checkbox.checked) {
        content.classList.remove('disabled');
        content.querySelectorAll('input, select').forEach(input => {
            input.disabled = false;
        });
    } else {
        content.classList.add('disabled');
        content.querySelectorAll('input, select').forEach(input => {
            input.disabled = true;
        });
    }
}

// Handle Start Simulation
async function handleStartSimulation(e) {
    e.preventDefault();

    const modelName = document.getElementById('model_name').value;
    const statusDiv = document.getElementById('setupStatus');

    // Clear any previous status messages
    statusDiv.textContent = '';
    statusDiv.className = 'status-message';

    // Collect all form parameters dynamically
    const formParams = collectFormParameters();

    // Add required metadata
    const params = {
        username: state.currentUser,
        base_model: document.getElementById('base_model').value,
        model_name: modelName,
        map_name: document.getElementById('map_name')?.value || 'sandbox',
        ...formParams
    };

    console.log('Starting simulation with params:', params);

    try {
        const response = await fetch(`${API_BASE}/simulation/start`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params)
        });

        console.log('Start simulation response status:', response.status);
        const data = await response.json();
        console.log('Start simulation response data:', data);

        if (response.ok) {
            // Clear any previous job state
            stopProgressPolling();

            // Reset progress bar immediately
            document.getElementById('progressBar').style.width = '0%';
            document.getElementById('progressPercent').textContent = '0%';

            // Reset plot state and hide plots display
            plotState.allData = {};
            plotState.lastRowCount = 0;
            plotState.lastModified = null;
            plotState.availableMetrics = [];
            if (plotState.updateInterval) {
                clearInterval(plotState.updateInterval);
                plotState.updateInterval = null;
            }
            // Destroy existing charts to prevent memory leaks
            plotState.charts.forEach(chart => {
                if (chart && typeof chart.destroy === 'function') {
                    chart.destroy();
                }
            });
            plotState.charts = [];
            document.getElementById('noPlots').style.display = 'block';
            document.getElementById('plotsDisplay').style.display = 'none';

            // Reset maps display
            document.getElementById('noMaps').style.display = 'block';
            document.getElementById('mapsDisplay').style.display = 'none';
            // Clear map image sources to avoid showing stale images
            const geneticMapAnim = document.getElementById('geneticMapAnim');
            const genealogicalMapAnim = document.getElementById('genealogicalMapAnim');
            if (geneticMapAnim) geneticMapAnim.src = '';
            if (genealogicalMapAnim) genealogicalMapAnim.src = '';

            // Set new job
            state.currentJobId = data.id;
            state.currentModel = modelName;

            console.log('Job created:', { jobId: data.id, modelName: modelName });

            // Switch to progress tab immediately (don't show success message on setup page)
            switchTab('progress');
            startProgressPolling();
        } else {
            throw new Error(data.error || 'Failed to start simulation');
        }
    } catch (error) {
        console.error('Start simulation error:', error);
        // Provide helpful error message for common issues
        let errorMessage = error.message;
        if (error.message === 'Failed to fetch') {
            errorMessage = 'Cannot connect to server. Please ensure the server is running (go run drift.go -web)';
        }
        statusDiv.textContent = 'Error: ' + errorMessage;
        statusDiv.className = 'status-message error';
    }
}

// Save current model parameters
async function handleSaveModel(e) {
    e.preventDefault();

    const modelName = document.getElementById('model_name').value.trim();

    if (!state.currentUser) {
        alert('No user logged in');
        return;
    }

    if (!modelName) {
        alert('Please enter a model name');
        return;
    }

    const statusDiv = document.getElementById('setupStatus');

    try {
        // Collect all parameters from the form
        const parameters = collectFormParameters();

        console.log('Parameters to save:', parameters);

        const payload = {
            username: state.currentUser,
            model_name: modelName,
            parameters: parameters
        };

        console.log('Payload:', payload);
        console.log('Payload JSON:', JSON.stringify(payload));

        const response = await fetch(`${API_BASE}/models/save`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });

        const data = await response.json();
        if (data.success) {
            state.currentModel = modelName;  // Update state to reflect saved model
            statusDiv.textContent = `Model "${modelName}" saved successfully!`;
            statusDiv.className = 'status-message success';
            setTimeout(() => { statusDiv.textContent = ''; }, 3000);
        } else {
            throw new Error(data.message || 'Failed to save model');
        }
    } catch (error) {
        console.error('Save error:', error);
        statusDiv.textContent = 'Error saving model: ' + error.message;
        statusDiv.className = 'status-message error';
    }
}

// Reset form to base model defaults
async function handleResetToDefaults(e) {
    e.preventDefault();

    const baseModelId = document.getElementById('base_model').value;

    if (!baseModelId) {
        alert('Please select a base model first');
        return;
    }

    const statusDiv = document.getElementById('setupStatus');

    try {
        const response = await fetch(`${API_BASE}/basemodels/parameters?id=${baseModelId}`);
        if (response.ok) {
            const data = await response.json();
            if (data.success && data.parameters) {
                populateFormFromParameters(data.parameters);
                document.getElementById('model_name').value = baseModelId;
                statusDiv.textContent = `Reset to ${baseModelId} defaults`;
                statusDiv.className = 'status-message success';
                setTimeout(() => { statusDiv.textContent = ''; }, 3000);
            } else {
                throw new Error('Failed to load base model parameters');
            }
        } else {
            throw new Error('Base model not found');
        }
    } catch (error) {
        statusDiv.textContent = 'Error resetting: ' + error.message;
        statusDiv.className = 'status-message error';
    }
}

// Save current configuration as a new model
async function handleSaveAsNewModel(e) {
    e.preventDefault();

    if (!state.currentUser) {
        alert('No user logged in');
        return;
    }

    const newModelName = prompt('Enter name for new model:');
    if (!newModelName) return;

    const statusDiv = document.getElementById('setupStatus');

    try {
        // Collect all parameters from the form
        const parameters = collectFormParameters();

        // For now, we'll use the save endpoint to create a new model
        // This assumes the backend will create the model if it doesn't exist
        const response = await fetch(`${API_BASE}/models/save`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                username: state.currentUser,
                model_name: newModelName,
                parameters: parameters
            })
        });

        const data = await response.json();
        if (data.success) {
            statusDiv.textContent = `New model "${newModelName}" created!`;
            statusDiv.className = 'status-message success';

            // Reload the models list to show the new model
            await loadAvailableModels();

            // Select the new model
            document.getElementById('base_model').value = newModelName;
            state.currentModel = newModelName;

            setTimeout(() => { statusDiv.textContent = ''; }, 3000);
        } else {
            throw new Error(data.message || 'Failed to create model');
        }
    } catch (error) {
        statusDiv.textContent = 'Error creating model: ' + error.message;
        statusDiv.className = 'status-message error';
    }
}

// Collect all parameters from the form
function collectFormParameters() {
    const parameters = {};

    const parseIntSafe = (value, defaultValue) => {
        const parsed = parseInt(value);
        return isNaN(parsed) ? defaultValue : parsed;
    };

    const parseFloatSafe = (value, defaultValue) => {
        const parsed = parseFloat(value);
        return isNaN(parsed) ? defaultValue : parsed;
    };

    // Collect all input values from the setup form
    const form = document.getElementById('setupForm');
    if (!form) return parameters;

    // Get all inputs and selects from the form
    const inputs = form.querySelectorAll('input, select');

    inputs.forEach(input => {
        const id = input.id;
        if (!id) return;

        // Skip non-parameter inputs (handled separately)
        if (id === 'base_model' || id === 'model_name' || id === 'map_name') return;

        const metadata = PARAMETER_METADATA[id];
        const fieldType = metadata?.type || determineFieldType(id, input.value, metadata);

        if (input.type === 'checkbox') {
            parameters[id] = input.checked ? 1 : 0;
        } else if (input.tagName === 'SELECT') {
            // Selects with dropdown type that have numeric values
            if (fieldType === 'dropdown' && metadata?.options) {
                parameters[id] = parseIntSafe(input.value, 0);
            } else {
                // Keep as string if not a numeric dropdown
                parameters[id] = input.value;
            }
        } else if (fieldType === 'float') {
            parameters[id] = parseFloatSafe(input.value, 0);
        } else {
            parameters[id] = parseIntSafe(input.value, 0);
        }
    });

    return parameters;
}

// Check for Existing Simulation
async function checkForExistingSimulation() {
    try {
        const response = await fetch(`${API_BASE}/simulation/status?username=${state.currentUser}`);
        const data = await response.json();

        if (data.jobs && data.jobs.length > 0) {
            // Find the first running or queued job
            const activeJob = data.jobs.find(job => job.status === 'running' || job.status === 'queued');

            if (activeJob) {
                state.currentJobId = activeJob.id;
                state.currentModel = activeJob.model_name;
                startProgressPolling();
            } else {
                // Use the most recent job for display (even if completed)
                const latestJob = data.jobs[0];
                state.currentJobId = latestJob.id;
                state.currentModel = latestJob.model_name;

                // Don't poll if it's completed
                if (latestJob.status === 'completed') {
                    checkProgress(); // Check once to update display
                }
            }
        }
    } catch (error) {
        console.error('Failed to check for existing simulation:', error);
    }
}

// Check Progress
async function checkProgress() {
    console.log('checkProgress called with:', { user: state.currentUser, model: state.currentModel, jobId: state.currentJobId });

    if (!state.currentUser || !state.currentModel) {
        console.log('checkProgress: Missing user or model');
        return;
    }

    try {
        // First check job status to get current state
        const statusResponse = await fetch(`${API_BASE}/simulation/status?username=${state.currentUser}`);
        let jobStatus = null;

        console.log('Status response status:', statusResponse.status);
        if (statusResponse.ok) {
            const statusData = await statusResponse.json();
            console.log('Status data:', statusData);
            if (statusData.jobs && statusData.jobs.length > 0) {
                // Find the current job
                jobStatus = statusData.jobs.find(job => job.id === state.currentJobId);
                console.log('Found job status:', jobStatus);
            }
        }

        // Then check progress.json
        const response = await fetch(`${API_BASE}/progress/${state.currentUser}/${state.currentModel}`);
        console.log('Progress response status:', response.status);

        if (response.ok) {
            const progress = await response.json();
            console.log('Progress data:', progress);
            updateProgressDisplay(progress, jobStatus);
        } else {
            console.log('No progress file yet. Job status:', jobStatus);
            // No progress file yet - show job status if available
            if (jobStatus && (jobStatus.status === 'running' || jobStatus.status === 'queued')) {
                console.log('Showing initial job status');
                updateProgressDisplay({ status: jobStatus.status, current_generation: 0, total_generations: 0, population_size: 0 }, jobStatus);
            } else {
                console.log('No active job found, showing no simulation message');
                document.getElementById('noSimulation').style.display = 'block';
                document.getElementById('simulationProgress').style.display = 'none';
            }
        }
    } catch (error) {
        console.error('Failed to fetch progress:', error);
    }
}

// Update Progress Display
function updateProgressDisplay(progress, jobStatus) {
    document.getElementById('noSimulation').style.display = 'none';
    document.getElementById('simulationProgress').style.display = 'block';

    // Use job status if available, otherwise use progress.json status
    const actualStatus = jobStatus ? jobStatus.status : progress.status;

    // Update stats
    document.getElementById('simStatus').textContent = actualStatus;
    document.getElementById('simStatus').className = `status-badge ${actualStatus}`;
    document.getElementById('currentGen').textContent = progress.current_generation || 0;
    document.getElementById('totalGen').textContent = progress.total_generations || 0;
    document.getElementById('popSize').textContent = progress.population_size || 0;

    // Show error message if failed
    const errorDiv = document.getElementById('errorMessage');
    if (actualStatus === 'failed') {
        console.log('Job failed. jobStatus:', jobStatus);
        if (jobStatus && jobStatus.error) {
            errorDiv.textContent = 'Error: ' + jobStatus.error;
            errorDiv.style.display = 'block';
        } else {
            errorDiv.textContent = 'Simulation failed (no error details available)';
            errorDiv.style.display = 'block';
        }
    } else {
        errorDiv.style.display = 'none';
    }

    // Update progress bar
    let percent = 0;
    if (progress.total_generations > 0) {
        percent = Math.round((progress.current_generation / progress.total_generations) * 100);
    }
    document.getElementById('progressBar').style.width = percent + '%';
    document.getElementById('progressPercent').textContent = percent + '%';

    // If completed or cancelled, stop polling and load results
    if (actualStatus === 'completed' || actualStatus === 'cancelled' || actualStatus === 'failed') {
        stopProgressPolling();
        if (actualStatus === 'completed') {
            // Results are available - both plots and maps can be loaded when user navigates to those tabs
        }
    }
}

// Start Progress Polling
function startProgressPolling() {
    // Check immediately
    checkProgress();

    // Then poll every 2 seconds
    if (state.progressInterval) {
        clearInterval(state.progressInterval);
    }
    state.progressInterval = setInterval(checkProgress, 2000);
}

// Stop Progress Polling
function stopProgressPolling() {
    if (state.progressInterval) {
        clearInterval(state.progressInterval);
        state.progressInterval = null;
    }
}

// Handle Cancel Simulation
async function handleCancelSimulation() {
    if (!state.currentJobId) return;

    if (!confirm('Are you sure you want to cancel this simulation?')) return;

    try {
        const response = await fetch(`${API_BASE}/simulation/cancel`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ job_id: state.currentJobId })
        });

        if (response.ok) {
            stopProgressPolling();

            // Update the display to show cancelled status immediately
            checkProgress();

            // Show success message
            const statusDiv = document.getElementById('setupStatus');
            statusDiv.textContent = 'Simulation cancelled successfully. You can start a new simulation.';
            statusDiv.className = 'status-message success';

            // Clear message after 5 seconds
            setTimeout(() => {
                statusDiv.textContent = '';
                statusDiv.className = 'status-message';
            }, 5000);
        } else {
            const data = await response.json();
            alert('Failed to cancel: ' + (data.error || 'Unknown error'));
        }
    } catch (error) {
        alert('Error cancelling simulation: ' + error.message);
    }
}

// Load Plots (CSV data - future: generate plot images)
async function loadPlots() {
    if (!state.currentUser || !state.currentModel) {
        console.log('loadPlots: Missing user or model');
        return;
    }

    console.log('loadPlots: Checking for results CSV');

    try {
        const response = await fetch(`${API_BASE}/results/${state.currentUser}/${state.currentModel}`);

        if (response.ok) {
            const data = await response.json();
            console.log('loadPlots: Files found:', data.files);

            // Check if results CSV exists
            const hasCSV = data.files && data.files.some(f => f.name.endsWith('_results.csv'));

            if (hasCSV) {
                document.getElementById('noPlots').style.display = 'none';
                document.getElementById('plotsDisplay').style.display = 'block';
                console.log('loadPlots: Results CSV available');
            } else {
                console.log('loadPlots: No results CSV found');
                document.getElementById('noPlots').style.display = 'block';
                document.getElementById('plotsDisplay').style.display = 'none';
            }
        } else {
            console.log('loadPlots: Failed to fetch results');
            document.getElementById('noPlots').style.display = 'block';
            document.getElementById('plotsDisplay').style.display = 'none';
        }
    } catch (error) {
        console.error('loadPlots: Error:', error);
        document.getElementById('noPlots').style.display = 'block';
        document.getElementById('plotsDisplay').style.display = 'none';
    }
}

// Load Maps (GIF animations)
async function loadMaps() {
    if (!state.currentUser || !state.currentModel) {
        console.log('loadMaps: Missing user or model');
        return;
    }

    console.log('loadMaps: Loading map animations for', state.currentUser, state.currentModel);

    try {
        const response = await fetch(`${API_BASE}/results/${state.currentUser}/${state.currentModel}`);

        if (response.ok) {
            const data = await response.json();
            console.log('loadMaps: Files found:', data.files);

            // Check specifically for map animation files
            const hasGeneticGif = data.files && data.files.some(f => f.name === 'genetic.gif');
            const hasGenealogicalGif = data.files && data.files.some(f => f.name === 'genealogical.gif');

            console.log('loadMaps: Has genetic.gif:', hasGeneticGif, 'Has genealogical.gif:', hasGenealogicalGif);

            if (hasGeneticGif || hasGenealogicalGif) {
                document.getElementById('noMaps').style.display = 'none';
                document.getElementById('mapsDisplay').style.display = 'block';

                // Load map animation images
                const baseUrl = `${API_BASE}/results/${state.currentUser}/${state.currentModel}`;

                if (hasGeneticGif) {
                    document.getElementById('geneticMapAnim').src = `${baseUrl}/genetic.gif?t=${Date.now()}`;
                    console.log('loadMaps: Loading genetic.gif');
                } else {
                    document.getElementById('geneticMapAnim').alt = 'Genetic map animation not generated';
                }

                if (hasGenealogicalGif) {
                    document.getElementById('genealogicalMapAnim').src = `${baseUrl}/genealogical.gif?t=${Date.now()}`;
                    console.log('loadMaps: Loading genealogical.gif');
                } else {
                    document.getElementById('genealogicalMapAnim').alt = 'Genealogical map animation not generated';
                }
            } else {
                console.log('loadMaps: No map animation files found');
                if (data.files && data.files.length > 0) {
                    console.log('loadMaps: Simulation completed but map animations not generated. Track on Map may not have been enabled.');
                }
                document.getElementById('noMaps').style.display = 'block';
                document.getElementById('mapsDisplay').style.display = 'none';
            }
        } else {
            console.log('loadMaps: Failed to fetch results');
            document.getElementById('noMaps').style.display = 'block';
            document.getElementById('mapsDisplay').style.display = 'none';
        }
    } catch (error) {
        console.error('loadMaps: Error loading maps:', error);
        document.getElementById('noMaps').style.display = 'block';
        document.getElementById('mapsDisplay').style.display = 'none';
    }
}

// Download File
function downloadFile(filename) {
    const url = `${API_BASE}/results/${state.currentUser}/${state.currentModel}/${filename}`;
    window.open(url, '_blank');
}

// Utility Functions
function showError(message) {
    alert(message);
}

function showSuccess(message) {
    alert(message);
}

// ============================================================================
// PLOTTING FUNCTIONALITY
// ============================================================================

// Plot state management
const plotState = {
    charts: [],
    lastRowCount: 0,
    lastModified: null,
    updateInterval: null,
    isPaused: false,
    allData: {}, // Stores all data keyed by metric name
    availableMetrics: []  // Populated dynamically from results CSV headers
};

// Metric definitions - maps CSV column names to display names
// These are dynamically discovered from results CSV headers
const METRIC_DEFINITIONS = {
    // Population metrics
    'n': { name: 'Population Size', group: 'Population', color: '#1f77b4' },
    'marrs': { name: 'Marriages', group: 'Population', color: '#ff7f0e' },
    'births': { name: 'Births', group: 'Population', color: '#2ca02c' },
    'randDs': { name: 'Random Deaths', group: 'Population', color: '#d62728' },
    'cullDs': { name: 'Culled Deaths', group: 'Population', color: '#9467bd' },

    // Ancestry metrics
    'YDes': { name: 'Y Descendants', group: 'Ancestry', color: '#8c564b' },
    'MtDes': { name: 'mt Descendants', group: 'Ancestry', color: '#e377c2' },
    'GeneaDes': { name: 'Genealogical Descendants', group: 'Ancestry', color: '#7f7f7f' },
    'GenetDes': { name: 'Genetic Descendants', group: 'Ancestry', color: '#bcbd22' },

    // DNA/Genetic metrics
    'nCents': { name: 'Centromeres', group: 'DNA', color: '#c5b0d5' },
    'nAlleles': { name: 'Alleles', group: 'DNA', color: '#c49c94' },
    'nBlocks': { name: 'Blocks', group: 'DNA', color: '#f7b6d2' },
    'TotFitness': { name: 'Total Fitness', group: 'DNA', color: '#9edae5' },
    'nMuts': { name: 'Number of Mutations', group: 'DNA', color: '#ff9896' },
    'PercSeedGenoRet': { name: 'Seed Genome %', group: 'DNA', color: '#dbdb8d' },
    'AvSeedGenoCov': { name: 'Seed Coverage', group: 'DNA', color: '#17becf' },
    'AvHet': { name: 'Average Heterozygosity', group: 'DNA', color: '#aec7e8' },
    'HomMin': { name: 'Homozygosity Min', group: 'DNA', color: '#ffbb78' },
    'HomMaj': { name: 'Homozygosity Maj', group: 'DNA', color: '#98df8a' },

    // Drift metrics
    'numSegSites': { name: 'Segregating Sites', group: 'Drift', color: '#1f77b4' },
    'allelesLost': { name: 'Alleles Lost', group: 'Drift', color: '#d62728' },
    'allelesFixed': { name: 'Alleles Fixed', group: 'Drift', color: '#2ca02c' },
    'avMinorFreq': { name: 'Avg Minor Allele Freq', group: 'Drift', color: '#ff7f0e' },
    'freqVariance': { name: 'Frequency Variance', group: 'Drift', color: '#9467bd' },
    'TajimaPi': { name: "Tajima's Pi", group: 'Drift', color: '#8c564b' },
    'WattersonsTheta': { name: "Watterson's Theta", group: 'Drift', color: '#e377c2' },
    'Fis': { name: 'Fis', group: 'Drift', color: '#7f7f7f' },
    'Ne': { name: 'Effective Pop Size (Ne)', group: 'Drift', color: '#bcbd22' },

    // SFS metrics
    'numSegSitesSFS': { name: 'SFS Seg Sites', group: 'SFS', color: '#1f77b4' },
    'singletons': { name: 'Singletons', group: 'SFS', color: '#ff7f0e' },
    'doubletons': { name: 'Doubletons', group: 'SFS', color: '#2ca02c' },
    'TajimasD': { name: "Tajima's D", group: 'SFS', color: '#d62728' },
    'FuLiD': { name: "Fu & Li's D", group: 'SFS', color: '#9467bd' },
    'FayWuH': { name: "Fay & Wu's H", group: 'SFS', color: '#8c564b' },
    'sfsMean': { name: 'SFS Mean', group: 'SFS', color: '#e377c2' },
    'sfsVar': { name: 'SFS Variance', group: 'SFS', color: '#7f7f7f' },
    'sfsSkew': { name: 'SFS Skew', group: 'SFS', color: '#bcbd22' },
    'sfsExtremes': { name: 'SFS Extremes', group: 'SFS', color: '#17becf' },

    // LD metrics
    'LD_short': { name: 'LD Short Range', group: 'LD', color: '#1f77b4' },
    'LD_medium': { name: 'LD Medium Range', group: 'LD', color: '#ff7f0e' },
    'LD_long': { name: 'LD Long Range', group: 'LD', color: '#2ca02c' },
    'HighLDPairs': { name: 'High LD Pairs', group: 'LD', color: '#d62728' },
    'LDDecaySlope': { name: 'LD Decay Slope', group: 'LD', color: '#9467bd' },

    // Coalescence metrics
    'meanCoalTime': { name: 'Mean Coalescence Time', group: 'Coalescence', color: '#1f77b4' },
    'maxCoalTime': { name: 'Max Coalescence Time', group: 'Coalescence', color: '#ff7f0e' },
    'coalPairs': { name: 'Coalescence Pairs', group: 'Coalescence', color: '#2ca02c' },
    'foundMRCA': { name: 'Found MRCA', group: 'Coalescence', color: '#d62728' },
    'tmrcaMean': { name: 'TMRCA Mean', group: 'Coalescence', color: '#9467bd' },
    'tmrcaMax': { name: 'TMRCA Max', group: 'Coalescence', color: '#8c564b' },
    'tmrcaVariance': { name: 'TMRCA Variance', group: 'Coalescence', color: '#e377c2' },
    'estimatedNe': { name: 'Estimated Ne', group: 'Coalescence', color: '#7f7f7f' },
    'tmrcaSampleSize': { name: 'TMRCA Sample Size', group: 'Coalescence', color: '#bcbd22' },

    // Y-Adam metrics
    'yadamID': { name: 'Y-Adam ID', group: 'Ancestry', color: '#17becf' },
    'yadamGensBack': { name: 'Y-Adam Gens Back', group: 'Ancestry', color: '#ff9896' },
    'yadamBirthYear': { name: 'Y-Adam Birth Year', group: 'Ancestry', color: '#c5b0d5' }
};

// Standard default chart - ALWAYS created for every run
const STANDARD_CHART = {
    id: 'standard-chart',
    title: 'Population Overview',
    metrics: ['n', 'births', 'randDs', 'cullDs'],  // CSV column names
    yAxisType: 'linear',
    isStandard: true
};

// Additional pre-defined charts - these are shown if the metrics exist in the data
const OPTIONAL_CHARTS = [
    {
        id: 'ancestry-chart',
        title: 'Ancestry Tracking',
        metrics: ['YDes', 'MtDes', 'GeneaDes', 'GenetDes'],
        yAxisType: 'linear'
    },
    {
        id: 'genetic-chart',
        title: 'Genetic Metrics',
        metrics: ['AvHet', 'nMuts', 'TotFitness'],
        yAxisType: 'linear'
    },
    {
        id: 'drift-chart',
        title: 'Drift Statistics',
        metrics: ['numSegSites', 'allelesLost', 'allelesFixed'],
        yAxisType: 'linear'
    },
    {
        id: 'sfs-chart',
        title: 'Site Frequency Spectrum',
        metrics: ['TajimasD', 'FuLiD', 'FayWuH'],
        yAxisType: 'linear'
    },
    {
        id: 'ld-chart',
        title: 'Linkage Disequilibrium',
        metrics: ['LD_short', 'LD_medium', 'LD_long'],
        yAxisType: 'linear'
    }
];

// Initialize plot controls
function initializePlotControls() {
    // Update mode selector
    document.getElementById('updateMode')?.addEventListener('change', (e) => {
        if (e.target.value === 'auto') {
            resumePlotUpdates();
        } else {
            pausePlotUpdates();
        }
    });

    // Pause button
    document.getElementById('pausePlots')?.addEventListener('click', togglePausePlots);

    // Refresh button
    document.getElementById('refreshPlots')?.addEventListener('click', () => {
        fetchPlotData(true); // Force refresh
    });

    // Add custom chart button
    document.getElementById('addCustomChart')?.addEventListener('click', showCustomChartDialog);
}

// Toggle pause/resume
function togglePausePlots() {
    if (plotState.isPaused) {
        resumePlotUpdates();
    } else {
        pausePlotUpdates();
    }
}

// Pause plot updates
function pausePlotUpdates() {
    plotState.isPaused = true;
    if (plotState.updateInterval) {
        clearInterval(plotState.updateInterval);
        plotState.updateInterval = null;
    }
    const btn = document.getElementById('pausePlots');
    if (btn) {
        btn.textContent = 'Resume';
        btn.classList.remove('btn-secondary');
        btn.classList.add('btn-primary');
    }
}

// Resume plot updates
function resumePlotUpdates() {
    plotState.isPaused = false;
    startPlotUpdatePolling();
    const btn = document.getElementById('pausePlots');
    if (btn) {
        btn.textContent = 'Pause';
        btn.classList.remove('btn-primary');
        btn.classList.add('btn-secondary');
    }
}

// Start polling for plot updates
function startPlotUpdatePolling() {
    // Clear any existing interval
    if (plotState.updateInterval) {
        clearInterval(plotState.updateInterval);
    }

    // Fetch immediately
    fetchPlotData();

    // Then poll every 3 seconds
    plotState.updateInterval = setInterval(() => {
        fetchPlotData();
    }, 3000);
}

// Stop polling for plot updates
function stopPlotUpdatePolling() {
    if (plotState.updateInterval) {
        clearInterval(plotState.updateInterval);
        plotState.updateInterval = null;
    }
}

// Fetch plot data from server
async function fetchPlotData(forceRefresh = false) {
    if (!state.currentUser || !state.currentModel) {
        return;
    }

    try {
        const url = `${API_BASE}/plot-data/${state.currentUser}/${state.currentModel}?from_row=${forceRefresh ? 0 : plotState.lastRowCount}`;

        const headers = {};
        if (!forceRefresh && plotState.lastModified) {
            headers['If-Modified-Since'] = plotState.lastModified;
        }

        const response = await fetch(url, { headers });

        if (response.status === 304) {
            // No changes
            return;
        }

        if (!response.ok) {
            if (response.status === 404) {
                // No data yet
                return;
            }
            throw new Error(`Failed to fetch plot data: ${response.statusText}`);
        }

        const data = await response.json();

        // If force refresh, reset data
        if (forceRefresh) {
            plotState.allData = {};
            plotState.lastRowCount = 0;
        }

        // Process new rows
        if (data.rows && data.rows.length > 0) {
            processPlotData(data.headers, data.rows);
            plotState.lastRowCount = data.total_rows;
            plotState.lastModified = data.modified;

            // Update all charts
            updateAllCharts();

            // Update last update time
            const now = new Date().toLocaleTimeString();
            const updateSpan = document.getElementById('lastPlotUpdate');
            if (updateSpan) {
                updateSpan.textContent = `Last update: ${now}`;
            }
        }
    } catch (error) {
        console.error('Error fetching plot data:', error);
    }
}

// Process plot data rows
function processPlotData(headers, rows) {
    // Initialize data arrays if first time
    if (Object.keys(plotState.allData).length === 0) {
        headers.forEach(header => {
            plotState.allData[header] = [];
        });
        // Store ALL metrics, not just ones with definitions
        plotState.availableMetrics = headers.filter(h => h !== 'run' && h !== 'year');

        console.log('Available CSV columns:', headers);
        console.log('Metrics with definitions:', headers.filter(h => METRIC_DEFINITIONS[h]));
        console.log('Available metrics for plotting:', plotState.availableMetrics);
    }

    // Append new data
    rows.forEach(row => {
        headers.forEach(header => {
            const value = row[header];
            if (value !== undefined && value !== '') {
                plotState.allData[header].push(parseFloat(value) || value);
            } else {
                plotState.allData[header].push(null);
            }
        });
    });

    console.log('Processed', rows.length, 'new rows. Total data points:', plotState.allData['year']?.length);
    console.log('Sample data - year:', plotState.allData['year']?.slice(0, 3));
    console.log('Sample data - n:', plotState.allData['n']?.slice(0, 3));
    console.log('Sample data - births:', plotState.allData['births']?.slice(0, 3));
}

// Load plots (called when switching to Plots tab)
async function loadPlots() {
    if (!state.currentUser || !state.currentModel) {
        console.log('loadPlots: Missing user or model');
        return;
    }

    console.log('loadPlots: Checking for results data');

    try {
        // Check if results file exists
        const response = await fetch(`${API_BASE}/results/${state.currentUser}/${state.currentModel}`);

        if (response.ok) {
            const data = await response.json();
            const hasCSV = data.files && data.files.some(f => f.name.endsWith('_results.csv'));

            if (hasCSV) {
                document.getElementById('noPlots').style.display = 'none';
                document.getElementById('plotsDisplay').style.display = 'block';

                // Initialize plot controls if not already done
                if (!plotState.updateInterval) {
                    initializePlotControls();
                }

                // Capture which metrics were selected in setup
                captureSelectedMetrics();

                // Create default charts if none exist
                if (plotState.charts.length === 0) {
                    createDefaultCharts();
                }

                // Start polling for updates
                if (!plotState.isPaused) {
                    startPlotUpdatePolling();
                }

                console.log('loadPlots: Results CSV available, charts initialized');
            } else {
                console.log('loadPlots: No results CSV found');
                document.getElementById('noPlots').style.display = 'block';
                document.getElementById('plotsDisplay').style.display = 'none';
            }
        } else {
            console.log('loadPlots: Failed to fetch results');
            document.getElementById('noPlots').style.display = 'block';
            document.getElementById('plotsDisplay').style.display = 'none';
        }
    } catch (error) {
        console.error('loadPlots: Error:', error);
        document.getElementById('noPlots').style.display = 'block';
        document.getElementById('plotsDisplay').style.display = 'none';
    }
}

// Capture which metrics were selected in the setup (no longer needed - charts are based on available data)
function captureSelectedMetrics() {
    // This function is kept for compatibility but no longer does anything
    // Charts are now dynamically created based on available metrics in the results CSV
}

// Create default charts based on selected metrics
function createDefaultCharts() {
    const container = document.getElementById('chartsContainer');
    if (!container) return;

    // Clear existing charts
    container.innerHTML = '';
    plotState.charts = [];

    // ALWAYS create the standard chart (N, births, deaths)
    createChart(STANDARD_CHART);

    // Create optional charts based on available metrics in the data
    OPTIONAL_CHARTS.forEach(chartConfig => {
        // Filter to only include metrics that exist in the available data
        const availableMetrics = chartConfig.metrics.filter(metric =>
            plotState.availableMetrics.includes(metric)
        );

        // Only create chart if at least one metric has non-zero data
        if (availableMetrics.length > 0) {
            createChart({
                ...chartConfig,
                metrics: availableMetrics
            });
        }
    });
}

// Create a single chart
function createChart(config) {
    const container = document.getElementById('chartsContainer');
    if (!container) return;

    // Create chart wrapper
    const chartWrapper = document.createElement('div');
    chartWrapper.className = 'chart-wrapper';
    chartWrapper.id = `wrapper-${config.id}`;

    // Create chart header
    const header = document.createElement('div');
    header.className = 'chart-header';

    const title = document.createElement('h3');
    title.textContent = config.title;
    header.appendChild(title);

    const controls = document.createElement('div');
    controls.className = 'chart-controls';

    // Add metric toggles
    config.metrics.forEach(metric => {
        const metricDef = METRIC_DEFINITIONS[metric];
        if (metricDef) {
            const label = document.createElement('label');
            label.className = 'metric-toggle';

            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.checked = true;
            checkbox.dataset.metric = metric;
            checkbox.dataset.chartId = config.id;
            checkbox.addEventListener('change', () => toggleMetricVisibility(config.id, metric, checkbox.checked));

            const span = document.createElement('span');
            span.textContent = metricDef.name;
            span.style.color = metricDef.color;

            label.appendChild(checkbox);
            label.appendChild(span);
            controls.appendChild(label);
        }
    });

    // Add settings button
    const settingsBtn = document.createElement('button');
    settingsBtn.className = 'btn-icon';
    settingsBtn.innerHTML = '⚙';
    settingsBtn.title = 'Chart Settings';
    settingsBtn.addEventListener('click', () => showChartSettings(config.id));
    controls.appendChild(settingsBtn);

    // Add delete button
    const deleteBtn = document.createElement('button');
    deleteBtn.className = 'btn-icon';
    deleteBtn.innerHTML = '×';
    deleteBtn.title = 'Remove Chart';
    deleteBtn.addEventListener('click', () => deleteChart(config.id));
    controls.appendChild(deleteBtn);

    header.appendChild(controls);
    chartWrapper.appendChild(header);

    // Create plot div
    const plotDiv = document.createElement('div');
    plotDiv.id = config.id;
    plotDiv.className = 'plot-container';
    chartWrapper.appendChild(plotDiv);

    container.appendChild(chartWrapper);

    // Initialize Plotly chart
    initializePlotlyChart(config);

    // Store chart configuration
    plotState.charts.push({
        id: config.id,
        config: config,
        visibleMetrics: new Set(config.metrics)
    });
}

// Initialize a Plotly chart
function initializePlotlyChart(config) {
    const traces = config.metrics.map(metric => {
        const metricDef = METRIC_DEFINITIONS[metric];
        return {
            x: [],
            y: [],
            name: metricDef ? metricDef.name : metric,
            type: 'scatter',
            mode: 'lines',
            line: { color: metricDef ? metricDef.color : undefined }
        };
    });

    const layout = {
        title: '',
        xaxis: { title: 'Year' },
        yaxis: {
            title: 'Value',
            type: config.yAxisType || 'linear'
        },
        margin: { t: 30, r: 30, b: 50, l: 60 },
        showlegend: true,
        legend: { orientation: 'h', y: -0.2 }
    };

    const plotConfig = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        displaylogo: false
    };

    Plotly.newPlot(config.id, traces, layout, plotConfig);
}

// Update all charts with current data
function updateAllCharts() {
    plotState.charts.forEach(chart => {
        updateChart(chart.id);
    });
}

// Update a specific chart
function updateChart(chartId) {
    const chart = plotState.charts.find(c => c.id === chartId);
    if (!chart) {
        console.log('Chart not found:', chartId);
        return;
    }

    // Use 'year' column for X-axis
    const xData = plotState.allData['year'] || plotState.allData['generation'] || plotState.allData['gen'] || [];

    console.log(`Updating chart ${chartId}:`, {
        metrics: chart.config.metrics,
        visibleMetrics: Array.from(chart.visibleMetrics),
        xDataLength: xData.length,
        sampleX: xData.slice(0, 3)
    });

    const updates = {
        x: [],
        y: []
    };

    chart.config.metrics.forEach((metric, index) => {
        if (chart.visibleMetrics.has(metric)) {
            const yData = plotState.allData[metric] || [];
            console.log(`  Metric ${metric}: ${yData.length} points, sample:`, yData.slice(0, 3));
            updates.x.push(xData);
            updates.y.push(yData);
        } else {
            updates.x.push([]);
            updates.y.push([]);
        }
    });

    console.log('Updating Plotly with:', updates);

    // Plotly.restyle requires trace indices - create array [0, 1, 2, ...] for all traces
    const traceIndices = chart.config.metrics.map((_, i) => i);
    Plotly.restyle(chartId, updates, traceIndices);
}

// Toggle metric visibility
function toggleMetricVisibility(chartId, metric, visible) {
    const chart = plotState.charts.find(c => c.id === chartId);
    if (!chart) return;

    if (visible) {
        chart.visibleMetrics.add(metric);
    } else {
        chart.visibleMetrics.delete(metric);
    }

    updateChart(chartId);
}

// Show chart settings dialog
function showChartSettings(chartId) {
    const chart = plotState.charts.find(c => c.id === chartId);
    if (!chart) return;

    const currentType = chart.config.yAxisType || 'linear';
    const newType = currentType === 'linear' ? 'log' : 'linear';

    // Update chart config
    chart.config.yAxisType = newType;

    // Update Plotly layout
    Plotly.relayout(chartId, {
        'yaxis.type': newType
    });
}

// Delete a chart
function deleteChart(chartId) {
    if (!confirm('Remove this chart?')) return;

    // Remove from DOM
    const wrapper = document.getElementById(`wrapper-${chartId}`);
    if (wrapper) {
        wrapper.remove();
    }

    // Remove from state
    const index = plotState.charts.findIndex(c => c.id === chartId);
    if (index !== -1) {
        plotState.charts.splice(index, 1);
    }
}

// Show custom chart creation dialog
function showCustomChartDialog() {
    if (plotState.availableMetrics.length === 0) {
        alert('No data available yet. Wait for simulation to generate data.');
        return;
    }

    // Create a simple dialog for selecting metrics
    const dialogHTML = `
        <div style="position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(0,0,0,0.5); z-index: 9999; display: flex; align-items: center; justify-content: center;" id="customChartDialog">
            <div style="background: white; padding: 30px; border-radius: 12px; max-width: 600px; max-height: 80vh; overflow-y: auto;">
                <h3 style="margin-top: 0; color: #667eea;">Create Custom Chart</h3>
                <div style="margin-bottom: 20px;">
                    <label style="display: block; margin-bottom: 5px; font-weight: 500;">Chart Name:</label>
                    <input type="text" id="customChartName" style="width: 100%; padding: 8px; border: 1px solid #ddd; border-radius: 5px;" placeholder="Enter chart name">
                </div>
                <div style="margin-bottom: 20px;">
                    <label style="display: block; margin-bottom: 10px; font-weight: 500;">Select Metrics (choose 1-6):</label>
                    <div id="metricCheckboxes" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px; max-height: 300px; overflow-y: auto; padding: 10px; border: 1px solid #ddd; border-radius: 5px;">
                        ${plotState.availableMetrics.map(metric => {
                            const def = METRIC_DEFINITIONS[metric];
                            const displayName = def ? def.name : metric;
                            return `
                                <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
                                    <input type="checkbox" value="${metric}" class="metric-checkbox">
                                    <span style="color: ${def ? def.color : '#333'}">${displayName}</span>
                                </label>
                            `;
                        }).join('')}
                    </div>
                </div>
                <div style="display: flex; gap: 10px; justify-content: flex-end;">
                    <button onclick="cancelCustomChart()" style="padding: 10px 20px; border: none; background: #6c757d; color: white; border-radius: 5px; cursor: pointer;">Cancel</button>
                    <button onclick="createCustomChart()" style="padding: 10px 20px; border: none; background: #667eea; color: white; border-radius: 5px; cursor: pointer;">Create Chart</button>
                </div>
            </div>
        </div>
    `;

    document.body.insertAdjacentHTML('beforeend', dialogHTML);
}

// Create custom chart from dialog
window.createCustomChart = function() {
    const chartName = document.getElementById('customChartName').value.trim();
    if (!chartName) {
        alert('Please enter a chart name');
        return;
    }

    const checkboxes = document.querySelectorAll('.metric-checkbox:checked');
    if (checkboxes.length === 0) {
        alert('Please select at least one metric');
        return;
    }
    if (checkboxes.length > 6) {
        alert('Please select no more than 6 metrics for readability');
        return;
    }

    const selectedMetrics = Array.from(checkboxes).map(cb => cb.value);

    const chartId = 'custom-' + Date.now();
    createChart({
        id: chartId,
        title: chartName,
        metrics: selectedMetrics,
        yAxisType: 'linear'
    });

    updateChart(chartId);

    // Close dialog
    cancelCustomChart();
};

// Cancel custom chart dialog
window.cancelCustomChart = function() {
    const dialog = document.getElementById('customChartDialog');
    if (dialog) {
        dialog.remove();
    }
};
