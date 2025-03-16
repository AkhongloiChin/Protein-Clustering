# syntax=docker/dockerfile:1
ARG PYTHON_VERSION=3.12.7
FROM python:${PYTHON_VERSION}-slim AS base

# Prevent .pyc files and enable unbuffered logs
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# Create a non-privileged user
ARG UID=10001
ARG GID=10001
RUN groupadd --gid "${GID}" appgroup && \
    useradd --uid "${UID}" --gid "${GID}" --no-create-home --shell /sbin/nologin appuser

# Fix matplotlib permission issue
RUN mkdir -p /tmp/matplotlib && chown -R appuser:appgroup /tmp/matplotlib
ENV MPLCONFIGDIR=/tmp/matplotlib

# Install dependencies
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

# Copy application source code
COPY . .

# Set the correct ownership after copying
RUN chown -R appuser:appgroup /app

# Copy entrypoint script and set executable permissions
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Switch to non-root user
USER appuser

# Expose the application's port
EXPOSE 8000

# Set entrypoint
ENTRYPOINT ["sh", "/entrypoint.sh"]

# Default command
CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
