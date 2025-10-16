# Autism Spectrum Disorder (ASD) Treatment Database

**Welcome to the ASD Treatment Database Github! Our vision is to enhance mental health outcomes for individuals with Autism Spectrum Disorder and their families by synthesizing psychiatric treatment knowledge for healthcare professionals, patients, and families. We aim to provide a comprehensive location for users of all kinds that allows for summarization of critical information of up-to-date medical papers to streamline your needs.**

[![Live Website](https://img.shields.io/badge/Live-Website-blue?style=for-the-badge)](https://star.cs.wwu.edu)
[![Environment Config](https://img.shields.io/badge/Multi-Environment-green?style=for-the-badge)](#environment-configuration)
[![Production Ready](https://img.shields.io/badge/Production-Ready-success?style=for-the-badge)](#production-deployment)

## **Project Overview**

This comprehensive database system provides:

- **Research Synthesis**: Automated extraction and analysis of ASD treatment literature from PubMed
- **Treatment Database**: Structured database of evidence-based interventions
- **Interactive Interface**: React-based web application for easy user access
- **AI-Powered Analysis**: MedBERT integration for advanced text processing

**Live System**: Generously hosted by Western Washington University | Database: Neon PostgreSQL

View our official docs [here](https://tipsycanoe.github.io/Autism-Spectrum-Disorder-ASD-Treatment-Database/)!

We also feature:

- **Advanced Search**: Filter treatments by age, symptoms, medications
- **Evidence Synthesis**: Automated analysis of treatment effectiveness
- **AI-Powered**: MedBERT integration for intelligent text processing
- **Responsive Design**: Works on desktop, tablet, and mobile
- **Auto-Updates**: Scheduled PubMed data refresh
- **Fast Performance**: Optimized queries and caching

## **Quick Start (Local Development)**

### Prerequisites

- **Python 3.8+** ([Download](https://www.python.org/downloads/))
- **Node.js 16+** ([Download](https://nodejs.org/))
- **Git** ([Download](https://git-scm.com/))

### Command Line Setup

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

Open [http://localhost:3000](http://localhost:3000) to view the application.

### Environment Commands

```bash
# Local development (default)
./start_all_servers.sh
```

### Configuration Files

- `config/local.env` - Local development settings

📚 **Detailed Guide**: See [ENVIRONMENT_GUIDE.md](ENVIRONMENT_GUIDE.md)

## **Production Deployment**

### Automated Production Deployment

```bash
# One-command production deployment
./deploy-production.sh
```

This handles:

- Environment setup and validation
- Frontend build process
- Service restarts and health checks
- Nginx configuration reload

For a more detailed guide and to see other deployment methods, please visit our dedicated deployment [page](https://tipsycanoe.github.io/Autism-Spectrum-Disorder-ASD-Treatment-Database/).

## **Architecture**

### System Components

- **Frontend**: React.js with Tailwind CSS
- **Backend APIs**:
  - Flask (Python) - Main application logic
  - Express (Node.js) - Job scheduling and automation
- **Database**: Neon PostgreSQL (cloud-hosted)
- **AI/ML**: MedBERT integration for text analysis
- **Deployment**: Nginx + Gunicorn + Systemd (production)

## **Database & APIs**

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

## **Documentation**

- **[Environment Setup Guide](mydocs/docs/getting_started/ENVIRONMENT_GUIDE.md)** - Comprehensive environment configuration
- **[Production Deployment](mydocs/docs/getting_started/PRODUCTION_DEPLOYMENT.md)** - Server deployment procedures
- **[File Overview](mydocs/docs/advanced/ENVIRONMENT_FILES_OVERVIEW.md)** - Complete file reference

## **Contributing**

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Set up local environment (`./start_all_servers.sh`)
4. Make your changes and test thoroughly
5. Commit your changes (`git commit -m 'Add AmazingFeature'`)
6. Push to the branch (`git push origin feature/AmazingFeature`)
7. Open a Pull Request

## 📄 **License**

This project is part of academic research at Western Washington University. Please contact the maintainers for usage permissions.

## **Acknowledgments**

- **Western Washington University** - Infrastructure and hosting support
- **Neon Database** - Cloud database hosting
- **PubMed/NCBI** - Research literature access
- **Hugging Face** - MedBERT model hosting
