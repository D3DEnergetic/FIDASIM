/**
 * API Client Service
 *
 * Handles all HTTP requests to the FIDASIM backend API
 */

import axios, { AxiosError, AxiosInstance } from 'axios';
import type {
  User,
  LoginRequest,
  LoginResponse,
  Job,
  JobSubmit,
  JobStats,
  FileInfo,
  StorageUsage,
} from '@/types';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000';

class ApiClient {
  private client: AxiosInstance;

  constructor() {
    this.client = axios.create({
      baseURL: API_URL,
      headers: {
        'Content-Type': 'application/json',
      },
    });

    // Request interceptor to add auth token
    this.client.interceptors.request.use(
      (config) => {
        const token = localStorage.getItem('token');
        if (token) {
          config.headers.Authorization = `Bearer ${token}`;
        }
        return config;
      },
      (error) => Promise.reject(error)
    );

    // Response interceptor for error handling
    this.client.interceptors.response.use(
      (response) => response,
      (error: AxiosError) => {
        if (error.response?.status === 401) {
          // Unauthorized - clear token and redirect to login
          localStorage.removeItem('token');
          localStorage.removeItem('user');
          window.location.href = '/login';
        }
        return Promise.reject(error);
      }
    );
  }

  // Authentication
  async login(credentials: LoginRequest): Promise<LoginResponse> {
    const response = await this.client.post<LoginResponse>('/api/auth/login', credentials);
    return response.data;
  }

  async logout(): Promise<void> {
    await this.client.post('/api/auth/logout');
  }

  async getCurrentUser(): Promise<User> {
    const response = await this.client.get<User>('/api/auth/me');
    return response.data;
  }

  // Jobs
  async submitJob(jobData: JobSubmit): Promise<Job> {
    const response = await this.client.post<Job>('/api/jobs/submit', jobData);
    return response.data;
  }

  async getJobs(status?: string, limit = 100, offset = 0): Promise<Job[]> {
    const params: any = { limit, offset };
    if (status) params.status_filter = status;

    const response = await this.client.get<Job[]>('/api/jobs/', { params });
    return response.data;
  }

  async getJobStats(): Promise<JobStats> {
    const response = await this.client.get<JobStats>('/api/jobs/stats');
    return response.data;
  }

  async getJob(jobId: string): Promise<Job> {
    const response = await this.client.get<Job>(`/api/jobs/${jobId}`);
    return response.data;
  }

  async cancelJob(jobId: string): Promise<{ message: string }> {
    const response = await this.client.delete<{ message: string }>(`/api/jobs/${jobId}/cancel`);
    return response.data;
  }

  async getJobLogs(jobId: string, tail = 10000): Promise<{ logs: string }> {
    const response = await this.client.get<{ logs: string }>(`/api/jobs/${jobId}/logs`, {
      params: { tail },
    });
    return response.data;
  }

  // Files
  async uploadFile(file: File, description?: string): Promise<FileInfo> {
    const formData = new FormData();
    formData.append('file', file);
    if (description) formData.append('description', description);

    const response = await this.client.post<FileInfo>('/api/files/upload', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    });
    return response.data;
  }

  async uploadMultipleFiles(files: File[]): Promise<{ files: FileInfo[] }> {
    const formData = new FormData();
    files.forEach((file) => formData.append('files', file));

    const response = await this.client.post<{ files: FileInfo[] }>(
      '/api/files/upload-multiple',
      formData,
      {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
      }
    );
    return response.data;
  }

  async listFiles(): Promise<{ files: FileInfo[] }> {
    const response = await this.client.get<{ files: FileInfo[] }>('/api/files/list');
    return response.data;
  }

  async downloadFile(filename: string): Promise<Blob> {
    const response = await this.client.get(`/api/files/download/${filename}`, {
      responseType: 'blob',
    });
    return response.data;
  }

  async getFileInfo(filename: string): Promise<FileInfo> {
    const response = await this.client.get<FileInfo>(`/api/files/info/${filename}`);
    return response.data;
  }

  async deleteFile(filename: string): Promise<{ message: string }> {
    const response = await this.client.delete<{ message: string }>(`/api/files/${filename}`);
    return response.data;
  }

  async getStorageUsage(): Promise<StorageUsage> {
    const response = await this.client.get<StorageUsage>('/api/files/storage-usage');
    return response.data;
  }

  // Users
  async getUsers(limit = 100, offset = 0): Promise<User[]> {
    const response = await this.client.get<User[]>('/api/users/', {
      params: { limit, offset },
    });
    return response.data;
  }

  async createUser(userData: {
    username: string;
    email?: string;
    password: string;
    role?: string;
    max_cores?: number;
    max_jobs?: number;
  }): Promise<User> {
    const response = await this.client.post<User>('/api/users/', userData);
    return response.data;
  }

  async updateUser(
    userId: string,
    userData: Partial<{
      email: string;
      password: string;
      role: string;
      is_active: boolean;
      max_cores: number;
      max_jobs: number;
    }>
  ): Promise<User> {
    const response = await this.client.patch<User>(`/api/users/${userId}`, userData);
    return response.data;
  }

  async deleteUser(userId: string): Promise<{ message: string }> {
    const response = await this.client.delete<{ message: string }>(`/api/users/${userId}`);
    return response.data;
  }

  async activateUser(userId: string): Promise<{ message: string }> {
    const response = await this.client.post<{ message: string }>(`/api/users/${userId}/activate`);
    return response.data;
  }
}

export const api = new ApiClient();
