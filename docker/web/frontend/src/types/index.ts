// TypeScript type definitions for FIDASIM Web

export interface User {
  id: string;
  username: string;
  email?: string;
  role: 'admin' | 'user' | 'viewer';
  is_active: boolean;
  max_cores: number;
  max_jobs: number;
  created_at: string;
  updated_at: string;
  last_login?: string;
}

export interface LoginRequest {
  username: string;
  password: string;
}

export interface LoginResponse {
  access_token: string;
  token_type: string;
  user: User;
}

export type JobStatus = 'pending' | 'running' | 'completed' | 'failed' | 'cancelled';

export interface Job {
  id: string;
  user_id: string;
  run_id: string;
  status: JobStatus;
  config: Record<string, any>;
  cores_assigned: number;
  priority: number;
  input_files: string[];
  output_files: string[];
  error_message?: string;
  created_at: string;
  started_at?: string;
  completed_at?: string;
  metadata?: Record<string, any>;
}

export interface JobConfig {
  input_file: string;
  device?: string;
  result_dir?: string;
  tables_file?: string;
  passive_only?: boolean;
  cores: number;
}

export interface JobSubmit {
  run_id: string;
  config: JobConfig;
  priority?: number;
}

export interface JobStats {
  total: number;
  pending: number;
  running: number;
  completed: number;
  failed: number;
  cancelled: number;
}

export interface FileInfo {
  file_id?: string;
  filename: string;
  unique_filename?: string;
  size: number;
  size_human?: string;
  upload_date?: string;
  created?: string;
  modified?: string;
  path?: string;
  description?: string;
}

export interface StorageUsage {
  total_size: number;
  total_size_human: string;
  file_count: number;
  max_size: number;
  max_size_human: string;
}

export interface WebSocketMessage {
  type: 'connected' | 'job_update' | 'job_log' | 'notification' | 'pong';
  job_id?: string;
  status?: JobStatus;
  log?: string;
  notification_type?: 'info' | 'warning' | 'error' | 'success';
  title?: string;
  body?: string;
  message?: string;
  timestamp?: string;
  details?: Record<string, any>;
  user_id?: string;
  username?: string;
}

export interface ApiError {
  detail: string;
}
