# GenomicsToolkit Docker Image
# Multi-stage build for optimized production image

# Build stage
FROM python:3.11-slim as builder

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies for building
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    pkg-config \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create and activate virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt

# Production stage
FROM python:3.11-slim as production

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH" \
    GENOMICS_TOOLKIT_CONFIG="/app/config" \
    GENOMICS_TOOLKIT_DATA="/app/data" \
    GENOMICS_TOOLKIT_OUTPUT="/app/output"

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    libopenblas0 \
    liblapack3 \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy virtual environment from builder stage
COPY --from=builder /opt/venv /opt/venv

# Create application user and directories
RUN groupadd -r genomics && useradd -r -g genomics genomics && \
    mkdir -p /app /app/config /app/data /app/output && \
    chown -R genomics:genomics /app

# Set working directory
WORKDIR /app

# Copy application code
COPY --chown=genomics:genomics . .

# Install the package in development mode
RUN pip install -e .

# Create entrypoint script
RUN cat > /app/entrypoint.sh << 'EOF'
#!/bin/bash
set -e

# Initialize configuration if not exists
if [ ! -f "$GENOMICS_TOOLKIT_CONFIG/config.yaml" ]; then
    echo "Initializing configuration..."
    genomics-toolkit config --list-all > /dev/null
fi

# Create output directory if it doesn't exist
mkdir -p "$GENOMICS_TOOLKIT_OUTPUT"

# Execute the command
exec "$@"
EOF

RUN chmod +x /app/entrypoint.sh && chown genomics:genomics /app/entrypoint.sh

# Switch to non-root user
USER genomics

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD genomics-toolkit --help > /dev/null || exit 1

# Set entrypoint and default command
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["genomics-toolkit", "--help"]

# Labels for metadata
LABEL maintainer="your.email@example.com" \
      version="1.0.0" \
      description="GenomicsToolkit - Comprehensive bioinformatics pipeline" \
      org.opencontainers.image.title="GenomicsToolkit" \
      org.opencontainers.image.description="Comprehensive bioinformatics pipeline for DNA/RNA sequence analysis and variant calling" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.authors="Your Name <your.email@example.com>" \
      org.opencontainers.image.source="https://github.com/Febo2788/genomics-toolkit" \
      org.opencontainers.image.licenses="MIT"