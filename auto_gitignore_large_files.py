#!/usr/bin/env python3
"""
Auto Gitignore Large Files

This script automatically finds files larger than 100MB in the current directory
and adds them to .gitignore if they're not already being ignored.

Usage:
    python auto_gitignore_large_files.py [--size-limit SIZE_IN_MB] [--dry-run]
    
Examples:
    python auto_gitignore_large_files.py                    # Default 100MB limit
    python auto_gitignore_large_files.py --size-limit 50   # 50MB limit
    python auto_gitignore_large_files.py --dry-run         # Show what would be added without modifying
"""

import os
import argparse
import subprocess
from pathlib import Path
import sys


def get_git_root():
    """Get the git root directory."""
    try:
        result = subprocess.run(['git', 'rev-parse', '--show-toplevel'], 
                              capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        print("Error: Not in a git repository")
        sys.exit(1)


def is_git_ignored(file_path, git_root):
    """Check if a file is already ignored by git."""
    try:
        # Get relative path from git root
        rel_path = os.path.relpath(file_path, git_root)
        result = subprocess.run(['git', 'check-ignore', rel_path], 
                              capture_output=True, text=True)
        return result.returncode == 0
    except Exception:
        return False


def find_large_files(directory, size_limit_mb):
    """Find files larger than the specified size limit."""
    size_limit_bytes = size_limit_mb * 1024 * 1024
    large_files = []
    
    for root, dirs, files in os.walk(directory):
        # Skip .git directory
        if '.git' in dirs:
            dirs.remove('.git')
        
        for file in files:
            file_path = os.path.join(root, file)
            try:
                file_size = os.path.getsize(file_path)
                if file_size > size_limit_bytes:
                    large_files.append((file_path, file_size))
            except (OSError, FileNotFoundError):
                continue
    
    return large_files


def read_gitignore(gitignore_path):
    """Read existing .gitignore file."""
    if not os.path.exists(gitignore_path):
        return []
    
    with open(gitignore_path, 'r') as f:
        return [line.strip() for line in f.readlines()]


def write_gitignore(gitignore_path, lines):
    """Write lines to .gitignore file."""
    with open(gitignore_path, 'w') as f:
        for line in lines:
            f.write(line + '\n')


def add_to_gitignore(gitignore_path, file_path, git_root):
    """Add a file path to .gitignore."""
    # Get relative path from git root
    rel_path = os.path.relpath(file_path, git_root)
    
    # Read existing .gitignore
    lines = read_gitignore(gitignore_path)
    
    # Check if already in .gitignore
    if rel_path in lines:
        return False
    
    # Add comment and path
    lines.append(f"# Large file: {rel_path}")
    lines.append(rel_path)
    
    # Write back to .gitignore
    write_gitignore(gitignore_path, lines)
    return True


def format_file_size(size_bytes):
    """Format file size in human readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def main():
    parser = argparse.ArgumentParser(description='Automatically add large files to .gitignore')
    parser.add_argument('--size-limit', type=int, default=100, 
                       help='Size limit in MB (default: 100)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be added without modifying .gitignore')
    
    args = parser.parse_args()
    
    # Get git root and .gitignore path
    git_root = get_git_root()
    gitignore_path = os.path.join(git_root, '.gitignore')
    
    print(f"Scanning for files larger than {args.size_limit}MB...")
    print(f"Git root: {git_root}")
    print(f".gitignore: {gitignore_path}")
    print("-" * 50)
    
    # Find large files
    large_files = find_large_files(git_root, args.size_limit)
    
    if not large_files:
        print("No files larger than the specified limit found.")
        return
    
    # Sort by size (largest first)
    large_files.sort(key=lambda x: x[1], reverse=True)
    
    print(f"Found {len(large_files)} files larger than {args.size_limit}MB:")
    print()
    
    files_to_add = []
    
    for file_path, file_size in large_files:
        rel_path = os.path.relpath(file_path, git_root)
        size_str = format_file_size(file_size)
        
        # Check if already ignored
        if is_git_ignored(file_path, git_root):
            status = "✓ Already ignored"
        else:
            status = "✗ Not ignored"
            files_to_add.append((file_path, file_size))
        
        print(f"{rel_path:<50} {size_str:>10} {status}")
    
    print()
    
    if not files_to_add:
        print("All large files are already being ignored by git.")
        return
    
    print(f"Files to add to .gitignore: {len(files_to_add)}")
    
    if args.dry_run:
        print("\nDRY RUN - No changes will be made to .gitignore")
        print("Files that would be added:")
        for file_path, file_size in files_to_add:
            rel_path = os.path.relpath(file_path, git_root)
            size_str = format_file_size(file_size)
            print(f"  {rel_path} ({size_str})")
        return
    
    # Ask for confirmation
    response = input(f"\nAdd {len(files_to_add)} files to .gitignore? (y/N): ")
    if response.lower() != 'y':
        print("Operation cancelled.")
        return
    
    # Add files to .gitignore
    added_count = 0
    for file_path, file_size in files_to_add:
        if add_to_gitignore(gitignore_path, file_path, git_root):
            added_count += 1
            rel_path = os.path.relpath(file_path, git_root)
            print(f"Added: {rel_path}")
    
    print(f"\nSuccessfully added {added_count} files to .gitignore")
    print("Don't forget to commit the updated .gitignore file!")


if __name__ == "__main__":
    main()
