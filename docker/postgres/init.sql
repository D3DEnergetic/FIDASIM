-- FIDASIM Database Schema Initialization
-- PostgreSQL initialization script for job tracking and user management

-- Create extensions
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS "pgcrypto";

-- Create enum types
CREATE TYPE job_status AS ENUM ('pending', 'running', 'completed', 'failed', 'cancelled');
CREATE TYPE user_role AS ENUM ('admin', 'user', 'viewer');
CREATE TYPE log_level AS ENUM ('debug', 'info', 'warning', 'error', 'critical');

-- Users table
CREATE TABLE IF NOT EXISTS users (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    username VARCHAR(100) UNIQUE NOT NULL,
    email VARCHAR(255) UNIQUE,
    password_hash VARCHAR(255) NOT NULL,
    role user_role DEFAULT 'user',
    is_active BOOLEAN DEFAULT true,
    max_cores INTEGER DEFAULT 32,
    max_jobs INTEGER DEFAULT 10,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    last_login TIMESTAMP WITH TIME ZONE,
    metadata JSONB DEFAULT '{}'::jsonb
);

-- Jobs table
CREATE TABLE IF NOT EXISTS jobs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES users(id) ON DELETE CASCADE,
    run_id VARCHAR(255) NOT NULL,
    status job_status DEFAULT 'pending',
    cores_assigned INTEGER DEFAULT 4,
    priority INTEGER DEFAULT 5,

    -- Input configuration
    input_files JSONB DEFAULT '{}'::jsonb,
    config JSONB NOT NULL,
    namelist_content TEXT,

    -- Execution details
    container_id VARCHAR(100),
    node_assigned VARCHAR(100),
    pid INTEGER,

    -- Timing
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    queued_at TIMESTAMP WITH TIME ZONE,
    started_at TIMESTAMP WITH TIME ZONE,
    completed_at TIMESTAMP WITH TIME ZONE,

    -- Results
    output_files JSONB DEFAULT '[]'::jsonb,
    error_message TEXT,
    exit_code INTEGER,

    -- Performance metrics
    cpu_time_seconds FLOAT,
    memory_peak_mb FLOAT,
    particles_simulated BIGINT,

    -- Indexes for common queries
    CONSTRAINT jobs_run_id_unique UNIQUE(run_id)
);

CREATE INDEX idx_jobs_user_id ON jobs(user_id);
CREATE INDEX idx_jobs_status ON jobs(status);
CREATE INDEX idx_jobs_created_at ON jobs(created_at DESC);
CREATE INDEX idx_jobs_run_id ON jobs(run_id);

-- Parameter scans table
CREATE TABLE IF NOT EXISTS parameter_scans (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES users(id) ON DELETE CASCADE,
    scan_name VARCHAR(255) NOT NULL,
    description TEXT,
    base_config JSONB NOT NULL,
    parameter_name VARCHAR(255) NOT NULL,
    parameter_values JSONB NOT NULL,  -- Array of values
    total_runs INTEGER NOT NULL,
    completed_runs INTEGER DEFAULT 0,
    failed_runs INTEGER DEFAULT 0,
    status job_status DEFAULT 'pending',
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    completed_at TIMESTAMP WITH TIME ZONE,
    metadata JSONB DEFAULT '{}'::jsonb
);

CREATE INDEX idx_scans_user_id ON parameter_scans(user_id);
CREATE INDEX idx_scans_status ON parameter_scans(status);

-- Link table for scan jobs
CREATE TABLE IF NOT EXISTS scan_jobs (
    scan_id UUID REFERENCES parameter_scans(id) ON DELETE CASCADE,
    job_id UUID REFERENCES jobs(id) ON DELETE CASCADE,
    parameter_value JSONB NOT NULL,
    parameter_index INTEGER NOT NULL,
    PRIMARY KEY (scan_id, job_id)
);

-- Job logs table
CREATE TABLE IF NOT EXISTS job_logs (
    id BIGSERIAL PRIMARY KEY,
    job_id UUID REFERENCES jobs(id) ON DELETE CASCADE,
    timestamp TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    level log_level DEFAULT 'info',
    message TEXT NOT NULL,
    details JSONB DEFAULT '{}'::jsonb
);

CREATE INDEX idx_job_logs_job_id ON job_logs(job_id);
CREATE INDEX idx_job_logs_timestamp ON job_logs(timestamp DESC);
CREATE INDEX idx_job_logs_level ON job_logs(level);

-- Performance metrics table (detailed tracking)
CREATE TABLE IF NOT EXISTS performance_metrics (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    job_id UUID REFERENCES jobs(id) ON DELETE CASCADE,
    timestamp TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    metric_name VARCHAR(100) NOT NULL,
    metric_value FLOAT NOT NULL,
    unit VARCHAR(50),
    metadata JSONB DEFAULT '{}'::jsonb
);

CREATE INDEX idx_metrics_job_id ON performance_metrics(job_id);
CREATE INDEX idx_metrics_name ON performance_metrics(metric_name);
CREATE INDEX idx_metrics_timestamp ON performance_metrics(timestamp DESC);

-- File metadata table
CREATE TABLE IF NOT EXISTS files (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    job_id UUID REFERENCES jobs(id) ON DELETE CASCADE,
    user_id UUID REFERENCES users(id) ON DELETE CASCADE,
    filename VARCHAR(255) NOT NULL,
    filepath VARCHAR(500) NOT NULL,
    file_type VARCHAR(50),
    file_size BIGINT,
    checksum VARCHAR(64),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    metadata JSONB DEFAULT '{}'::jsonb
);

CREATE INDEX idx_files_job_id ON files(job_id);
CREATE INDEX idx_files_user_id ON files(user_id);
CREATE INDEX idx_files_filename ON files(filename);

-- System configuration table
CREATE TABLE IF NOT EXISTS system_config (
    key VARCHAR(100) PRIMARY KEY,
    value JSONB NOT NULL,
    description TEXT,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_by UUID REFERENCES users(id)
);

-- Audit log table
CREATE TABLE IF NOT EXISTS audit_log (
    id BIGSERIAL PRIMARY KEY,
    user_id UUID REFERENCES users(id),
    action VARCHAR(100) NOT NULL,
    resource_type VARCHAR(50),
    resource_id VARCHAR(100),
    details JSONB DEFAULT '{}'::jsonb,
    ip_address INET,
    user_agent TEXT,
    timestamp TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_audit_user_id ON audit_log(user_id);
CREATE INDEX idx_audit_timestamp ON audit_log(timestamp DESC);
CREATE INDEX idx_audit_action ON audit_log(action);

-- Views for common queries

-- Active jobs view
CREATE OR REPLACE VIEW active_jobs AS
SELECT
    j.id,
    j.run_id,
    j.status,
    j.cores_assigned,
    u.username,
    j.created_at,
    j.started_at,
    EXTRACT(EPOCH FROM (CURRENT_TIMESTAMP - j.started_at)) as runtime_seconds
FROM jobs j
LEFT JOIN users u ON j.user_id = u.id
WHERE j.status IN ('pending', 'running')
ORDER BY j.created_at DESC;

-- Job statistics view
CREATE OR REPLACE VIEW job_statistics AS
SELECT
    u.username,
    COUNT(*) as total_jobs,
    COUNT(*) FILTER (WHERE j.status = 'completed') as completed_jobs,
    COUNT(*) FILTER (WHERE j.status = 'failed') as failed_jobs,
    COUNT(*) FILTER (WHERE j.status = 'running') as running_jobs,
    AVG(j.cpu_time_seconds) as avg_cpu_time,
    SUM(j.cores_assigned) FILTER (WHERE j.status = 'running') as current_cores_used
FROM users u
LEFT JOIN jobs j ON u.id = j.user_id
GROUP BY u.username;

-- Functions and triggers

-- Update timestamp trigger
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ language 'plpgsql';

CREATE TRIGGER update_users_updated_at BEFORE UPDATE ON users
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

-- Function to calculate optimal core allocation
CREATE OR REPLACE FUNCTION calculate_optimal_cores(
    total_cores INTEGER,
    num_jobs INTEGER,
    scaling_efficiency JSONB DEFAULT '{}'::jsonb
) RETURNS INTEGER AS $$
DECLARE
    cores_per_job INTEGER;
BEGIN
    -- Simple algorithm, can be enhanced with actual scaling data
    IF num_jobs = 1 THEN
        RETURN LEAST(32, total_cores);
    END IF;

    -- Basic distribution
    cores_per_job := total_cores / num_jobs;

    -- Cap at maximum effective cores (typically 32 for OpenMP)
    cores_per_job := LEAST(cores_per_job, 32);

    -- Ensure at least 1 core per job
    cores_per_job := GREATEST(cores_per_job, 1);

    RETURN cores_per_job;
END;
$$ LANGUAGE plpgsql;

-- Function to get queue position
CREATE OR REPLACE FUNCTION get_queue_position(job_uuid UUID)
RETURNS INTEGER AS $$
DECLARE
    position INTEGER;
BEGIN
    SELECT COUNT(*) + 1 INTO position
    FROM jobs
    WHERE status = 'pending'
    AND created_at < (SELECT created_at FROM jobs WHERE id = job_uuid);

    RETURN position;
END;
$$ LANGUAGE plpgsql;

-- Insert default admin user (password: admin123 - CHANGE THIS!)
INSERT INTO users (username, email, password_hash, role)
VALUES (
    'admin',
    'admin@fidasim.local',
    crypt('admin123', gen_salt('bf')),
    'admin'
) ON CONFLICT (username) DO NOTHING;

-- Insert default system configuration
INSERT INTO system_config (key, value, description)
VALUES
    ('max_total_cores', '128'::jsonb, 'Maximum total cores available in the system'),
    ('max_cores_per_job', '32'::jsonb, 'Maximum cores allowed per job'),
    ('default_cores_per_job', '4'::jsonb, 'Default cores assigned if not specified'),
    ('job_timeout_seconds', '86400'::jsonb, 'Maximum job runtime in seconds (24 hours)'),
    ('max_queued_jobs', '1000'::jsonb, 'Maximum number of jobs in queue'),
    ('enable_auto_scaling', 'true'::jsonb, 'Enable automatic core allocation optimization')
ON CONFLICT (key) DO NOTHING;

-- Permissions (for row-level security if needed)
ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
ALTER TABLE files ENABLE ROW LEVEL SECURITY;

-- Create policies for row-level security
CREATE POLICY jobs_user_policy ON jobs
    FOR ALL
    USING (user_id = current_setting('app.current_user_id', true)::uuid OR
           EXISTS (SELECT 1 FROM users WHERE id = current_setting('app.current_user_id', true)::uuid AND role = 'admin'));

CREATE POLICY files_user_policy ON files
    FOR ALL
    USING (user_id = current_setting('app.current_user_id', true)::uuid OR
           EXISTS (SELECT 1 FROM users WHERE id = current_setting('app.current_user_id', true)::uuid AND role = 'admin'));

-- Indexes for performance
CREATE INDEX idx_jobs_status_created ON jobs(status, created_at DESC);
CREATE INDEX idx_jobs_user_status ON jobs(user_id, status);
CREATE INDEX idx_logs_job_timestamp ON job_logs(job_id, timestamp DESC);

-- Grant permissions
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO fidasim;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO fidasim;
GRANT ALL PRIVILEGES ON ALL FUNCTIONS IN SCHEMA public TO fidasim;