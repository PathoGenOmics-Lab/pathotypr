// ============================================================================
// Configuration and Constants
// ============================================================================

// Step definitions for each command
export const processSteps = {
  train: ['Loading data', 'Extracting k-mers', 'Training model', 'Saving'],
  predict: ['Loading model', 'Processing sequences', 'Classifying'],
  classify: ['Loading markers', 'Building k-mers', 'Classifying genomes'],
  'split-fastq': ['Indexing reference', 'Processing reads', 'Genotyping'],
  match: ['Indexing references', 'Processing reads', 'Matching']
};

// Panel information for navigation
export const panelInfo = {
  'home': {
    iconPath: '<path d="M3 9l9-7 9 7v11a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"/><polyline points="9 22 9 12 15 12 15 22"/>',
    title: 'Home',
    description: 'Welcome to Pathotypr - Select a tool to get started'
  },
  'train': {
    iconPath: '<path d="M12 20V10"/><path d="M18 20V4"/><path d="M6 20v-4"/>',
    title: 'Train Model',
    description: 'Build a Random Forest classifier from labeled sequences'
  },
  'predict': {
    iconPath: '<circle cx="12" cy="12" r="10"/><path d="M12 16v-4"/><path d="M12 8h.01"/>',
    title: 'Predict',
    description: 'Classify sequences using a trained model'
  },
  'classify': {
    iconPath: '<path d="M9 3H5a2 2 0 0 0-2 2v4m6-6h10a2 2 0 0 1 2 2v4M9 3v18m0 0h10a2 2 0 0 0 2-2v-4M9 21H5a2 2 0 0 1-2-2v-4"/>',
    title: 'Classify',
    description: 'Genotype genomes based on SNP markers'
  },
  'split-fastq': {
    iconPath: '<path d="M14.5 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V7.5L14.5 2z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/>',
    title: 'Split FASTQ',
    description: 'Alignment-free genotyping from raw reads'
  },
  'match': {
    iconPath: '<circle cx="11" cy="11" r="8"/><path d="m21 21-4.35-4.35"/>',
    title: 'Match Reference',
    description: 'Find best matching reference for each sample'
  }
};

// Default form parameters
export const defaultParameters = {
  train: {
    'train-kmer': 21,
    'train-split': 0.2,
    'train-threads': '',
    'train-cv-folds': '',
    'train-max-depth': 20,
    'train-min-leaf': 5
  },
  predict: {
    'predict-threads': ''
  },
  classify: {
    'classify-kmer': 31,
    'classify-threads': '',
    'classify-nested': false
  },
  splitfq: {
    'splitfq-kmer': 21,
    'splitfq-min-depth': 10,
    'splitfq-min-alt': 95,
    'splitfq-threads': '',
    'splitfq-paired': false,
    'splitfq-nested': false
  },
  match: {
    'match-kmer': 21,
    'match-threads': '',
    'match-max-ref-freq': 0,
    'match-coarse-top-n': 128,
    'match-coarse-stride': 8,
    'match-early-stop-confidence': 0,
    'match-early-stop-min-kmers': 1000000,
    'match-strict-percentages': true,
    'match-auto-tune': true
  }
};

// Lineage-specific colors for visualizations
export const lineageColors = {
  'A1': '#d1ae00',
  'A2': '#8ef5c8',
  'A3': '#73c2ff',
  'A4': '#ff9cdb',
  'L1': '#ff3091',
  'L2': '#001aff',
  'L3': '#8a0bd2',
  'L4': '#ff0000',
  'L5': '#995200',
  'L6': '#1eb040',
  'L7': '#fbff00',
  'L8': '#ff9d00',
  'L9': '#37ff30',
  'L10': '#8fbda1'
};

// Fallback colors for charts
export const fallbackColors = [
  '#14b8a6', '#8b5cf6', '#f59e0b', '#ef4444', '#3b82f6',
  '#ec4899', '#10b981', '#6366f1', '#f97316', '#06b6d4',
  '#84cc16', '#a855f7', '#22c55e', '#e11d48', '#0ea5e9'
];

// M. tuberculosis genotyping data URLs (markers + reference genome)
export const demoGenotypingData = {
  markers: {
    url: 'https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1',
    filename: 'pathotypr_lineage_markers_v1.0.0.tsv'
  },
  reference: {
    url: 'https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1',
    filename: 'reference.fasta'
  }
};

export const lineageMarkerData = {
  markers: {
    url: 'https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1',
    filename: 'pathotypr_lineage_markers_v1.0.0.tsv'
  },
  reference: {
    url: 'https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1',
    filename: 'MTB_ancestor_reference.fasta'
  }
};

export const drMarkerData = {
  markers: {
    url: 'https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1',
    filename: 'pathotypr_dr_markers_v1.0.0.tsv'
  },
  reference: {
    url: 'https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1',
    filename: 'MTB_ancestor_reference.fasta'
  }
};

// Pre-trained MTBC model for Predict
export const mtbcModel = {
  url: 'https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1',
  filename: 'pathotypr_rf_model_v1.0.0.pathotypr'
};

// UI timing constants (ms)
export const TIMING = {
  SCROLL_SETTLE: 100,        // Wait for layout before scrollIntoView
  CHART_RESIZE: 100,         // Wait for container resize before recreating charts
  TOAST_DURATION: 4000,      // Success toast + confetti visible time
  ANIMATION_FEEDBACK: 2000,  // File-generated / extension-error highlight
  FIELD_MISSING: 3000,       // Missing field highlight duration
  FULLSCREEN_EXIT: 400,      // Fullscreen exit CSS transition
  DROP_FEEDBACK: 500,        // Drop-received animation
  PROGRESS_HIDE: 500,        // Progress bar hide after completion
  SYSTEM_USAGE_POLL: 1000,   // Poll backend app usage continuously
  FORM_SAVE_DEBOUNCE: 500,   // localStorage save debounce
  FILTER_DEBOUNCE: 150,      // Results filter debounce
};

// Panel to category mapping for accent colors
export const categoryMap = {
  'home': null,
  'train': 'ml',
  'predict': 'ml',
  'classify': 'genotyping',
  'split-fastq': 'genotyping',
  'match': 'utils'
};
