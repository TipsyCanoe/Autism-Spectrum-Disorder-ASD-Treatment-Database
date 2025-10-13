# Autism Spectrum Disorder (ASD) Treatment Database

**Enhancing mental health outcomes for individuals with Autism Spectrum Disorder and their families by synthesizing psychiatric treatment knowledge for healthcare professionals, patients, and families.**

[![Live Website](https://img.shields.io/badge/Live-Website-blue?style=for-the-badge)](https://star.cs.wwu.edu)
[![Environment Config](https://img.shields.io/badge/Multi-Environment-green?style=for-the-badge)](#environment-configuration)
[![Production Ready](https://img.shields.io/badge/Production-Ready-success?style=for-the-badge)](#production-deployment)

## 🎯 **Project Overview**

This comprehensive database system provides:
- **Research Synthesis**: Automated extraction and analysis of ASD treatment literature from PubMed
- **Treatment Database**: Structured database of evidence-based interventions
- **Interactive Interface**: React-based web application for healthcare professionals
- **AI-Powered Analysis**: MedBERT integration for advanced text processing

**Live System**: Hosted by Western Washington University | Database: Neon PostgreSQL

## 🚀 **Quick Start (Local Development)**

### Prerequisites

- **Python 3.8+** ([Download](https://www.python.org/downloads/))
- **Node.js 16+** ([Download](https://nodejs.org/))
- **Git** ([Download](https://git-scm.com/))

### One-Command Setup

```bash
# Clone and setup
git clone https://github.com/TipsyCanoe/Autism-Spectrum-Disorder-ASD-Treatment-Database.git
cd Autism-Spectrum-Disorder-ASD-Treatment-Database

# Install dependencies and start
python3 -m venv venv && source venv/bin/activate
pip install -r backend/requirements.txt
cd frontend/testing-website && npm install && cd ../..

# Start all services
./start_all_servers.sh
```

**That's it!** 🎉 Open [http://localhost:3000](http://localhost:3000) to view the application.

### Services Started
- **Frontend**: React app on port 3000
- **Python API**: Flask backend on port 5000  
- **Node.js API**: Job scheduler on port 5001

## 🔧 **Environment Configuration**

This project supports multiple environments with automatic configuration:

### Available Environments
- **`local`** (default): Development with ports 3000, 5000, 5001
- **`staging`**: Pre-production testing environment
- **`production`**: Live deployment environment

### Environment Commands
```bash
# Local development (default)
./start_all_servers.sh

# Staging environment
ENVIRONMENT=staging ./start_all_servers.sh

# Load specific environment manually
source load-env.sh  # Loads local by default
ENVIRONMENT=production source load-env.sh
```

### Configuration Files
- `config/local.env` - Local development settings
- `config/staging.env` - Staging environment template
- `config/production.env.template` - Production template (copy to `production.env`)

📚 **Detailed Guide**: See [ENVIRONMENT_GUIDE.md](ENVIRONMENT_GUIDE.md)

## 🏗️ **Production Deployment**

### Automated Production Deployment
```bash
# One-command production deployment
./deploy-production.sh
```

This handles:
- ✅ Environment setup and validation
- ✅ Frontend build process
- ✅ Service restarts and health checks
- ✅ Nginx configuration reload

📚 **Production Guide**: See [PRODUCTION_DEPLOYMENT.md](PRODUCTION_DEPLOYMENT.md)

## 🧪 **Testing**

**Run all tests:**
```bash
./run_all_tests.sh
```

**Individual test suites:**
```bash
# Frontend tests
cd frontend/testing-website && npm test

# Backend tests  
cd backend/tests && ./run_tests.sh
```

## 🏛️ **Architecture**

### System Components
- **Frontend**: React.js with Tailwind CSS
- **Backend APIs**: 
  - Flask (Python) - Main application logic
  - Express (Node.js) - Job scheduling and automation
- **Database**: Neon PostgreSQL (cloud-hosted)
- **AI/ML**: MedBERT integration for text analysis
- **Deployment**: Nginx + Gunicorn + Systemd (production)

### Data Pipeline
1. **PubMed API Integration** → Automated literature extraction
2. **MedBERT Processing** → AI-powered text analysis and classification
3. **Database Storage** → Structured treatment and outcome data
4. **Web Interface** → Healthcare professional access and search

## 📊 **Features**

- **🔍 Advanced Search**: Filter treatments by age, symptoms, medications
- **📈 Evidence Synthesis**: Automated analysis of treatment effectiveness  
- **🤖 AI-Powered**: MedBERT integration for intelligent text processing
- **📱 Responsive Design**: Works on desktop, tablet, and mobile
- **🔄 Auto-Updates**: Scheduled PubMed data refresh
- **⚡ Fast Performance**: Optimized queries and caching

## 🗄️ **Database & APIs**

### PubMed Integration
- **Automated Extraction**: `pubmed_API_data.py`, `pubmed_API_ASD_data.py`
- **Manual Updates**: Run `API_JOB.py` or use web interface
- **Scheduling**: Configure in `/backend/scheduler.js` ([cron reference](https://crontab.guru))

### Database
- **Host**: Neon PostgreSQL (cloud)
- **Content**: Treatment studies, outcomes, patient demographics
- **Updates**: Automated nightly refresh from PubMed

### MedBERT/LLM Usage
```bash
# Activate environment
source venv/bin/activate

# Install ML dependencies
pip install -r FineTunedLLM/requirements.txt

# Run analysis scripts
python FineTunedLLM/MedBERT.py
```

## 📚 **Documentation**

- **[Environment Setup Guide](ENVIRONMENT_GUIDE.md)** - Comprehensive environment configuration
- **[Production Deployment](PRODUCTION_DEPLOYMENT.md)** - Server deployment procedures  
- **[File Overview](ENVIRONMENT_FILES_OVERVIEW.md)** - Complete file reference
- **Code Documentation** - Inline comments throughout codebase

## 🤝 **Contributing**

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Set up local environment (`./start_all_servers.sh`)
4. Make your changes and test thoroughly
5. Commit your changes (`git commit -m 'Add AmazingFeature'`)
6. Push to the branch (`git push origin feature/AmazingFeature`)
7. Open a Pull Request

## 📄 **License**

This project is part of academic research at Western Washington University. Please contact the maintainers for usage permissions.

## 🏥 **Acknowledgments**

- **Western Washington University** - Infrastructure and hosting support
- **Neon Database** - Cloud database hosting
- **PubMed/NCBI** - Research literature access
- **Hugging Face** - MedBERT model hosting
