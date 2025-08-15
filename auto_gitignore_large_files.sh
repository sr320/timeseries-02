#!/bin/bash

# Auto Gitignore Large Files - Bash Version
# This script automatically finds files larger than 100MB and adds them to .gitignore

set -e

# Default values
SIZE_LIMIT_MB=100
DRY_RUN=false
GIT_ROOT=""
GITIGNORE_PATH=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    -s, --size-limit SIZE    Size limit in MB (default: 100)
    -d, --dry-run           Show what would be added without modifying .gitignore
    -h, --help              Show this help message

Examples:
    $0                        # Default 100MB limit
    $0 -s 50                 # 50MB limit
    $0 --dry-run             # Show what would be added without modifying

EOF
}

# Function to parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -s|--size-limit)
                SIZE_LIMIT_MB="$2"
                shift 2
                ;;
            -d|--dry-run)
                DRY_RUN=true
                shift
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
    done
}

# Function to get git root directory
get_git_root() {
    if ! GIT_ROOT=$(git rev-parse --show-toplevel 2>/dev/null); then
        print_error "Not in a git repository"
        exit 1
    fi
    GITIGNORE_PATH="$GIT_ROOT/.gitignore"
}

# Function to check if file is already ignored by git
is_git_ignored() {
    local file_path="$1"
    local rel_path
    
    rel_path=$(realpath --relative-to="$GIT_ROOT" "$file_path" 2>/dev/null || echo "$file_path")
    
    if git check-ignore "$rel_path" >/dev/null 2>&1; then
        return 0  # File is ignored
    else
        return 1  # File is not ignored
    fi
}

# Function to format file size
format_file_size() {
    local size_bytes="$1"
    local size_mb=$((size_bytes / 1024 / 1024))
    local size_gb=$((size_bytes / 1024 / 1024 / 1024))
    
    if [ $size_gb -gt 0 ]; then
        echo "${size_gb}GB"
    elif [ $size_mb -gt 0 ]; then
        echo "${size_mb}MB"
    else
        echo "${size_bytes}B"
    fi
}

# Function to find large files
find_large_files() {
    local size_limit_bytes=$((SIZE_LIMIT_MB * 1024 * 1024))
    local large_files=()
    
    print_info "Scanning for files larger than ${SIZE_LIMIT_MB}MB..."
    
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            local file_size
            file_size=$(stat -c%s "$file" 2>/dev/null || echo "0")
            
            if [ "$file_size" -gt "$size_limit_bytes" ]; then
                large_files+=("$file:$file_size")
            fi
        fi
    done < <(find "$GIT_ROOT" -type f -size +${SIZE_LIMIT_MB}M -not -path "$GIT_ROOT/.git/*" -print0 2>/dev/null)
    
    echo "${large_files[@]}"
}

# Function to add file to .gitignore
add_to_gitignore() {
    local file_path="$1"
    local rel_path
    
    rel_path=$(realpath --relative-to="$GIT_ROOT" "$file_path" 2>/dev/null || echo "$file_path")
    
    # Check if already in .gitignore
    if grep -q "^$rel_path$" "$GITIGNORE_PATH" 2>/dev/null; then
        return 1
    fi
    
    # Add to .gitignore
    echo "" >> "$GITIGNORE_PATH"
    echo "# Large file: $rel_path" >> "$GITIGNORE_PATH"
    echo "$rel_path" >> "$GITIGNORE_PATH"
    
    return 0
}

# Main function
main() {
    parse_args "$@"
    
    print_info "Auto Gitignore Large Files"
    print_info "Size limit: ${SIZE_LIMIT_MB}MB"
    print_info "Git root: $GIT_ROOT"
    print_info ".gitignore: $GITIGNORE_PATH"
    echo "----------------------------------------"
    
    # Get git root
    get_git_root
    
    # Find large files
    local large_files
    large_files=($(find_large_files))
    
    if [ ${#large_files[@]} -eq 0 ]; then
        print_success "No files larger than ${SIZE_LIMIT_MB}MB found."
        exit 0
    fi
    
    print_info "Found ${#large_files[@]} files larger than ${SIZE_LIMIT_MB}MB:"
    echo
    
    local files_to_add=()
    
    # Process each large file
    for file_info in "${large_files[@]}"; do
        IFS=':' read -r file_path file_size <<< "$file_info"
        local rel_path
        rel_path=$(realpath --relative-to="$GIT_ROOT" "$file_path" 2>/dev/null || echo "$file_path")
        local size_str
        size_str=$(format_file_size "$file_size")
        
        # Check if already ignored
        if is_git_ignored "$file_path"; then
            printf "%-50s %10s %s\n" "$rel_path" "$size_str" "✓ Already ignored"
        else
            printf "%-50s %10s %s\n" "$rel_path" "$size_str" "✗ Not ignored"
            files_to_add+=("$file_path:$file_size")
        fi
    done
    
    echo
    
    if [ ${#files_to_add[@]} -eq 0 ]; then
        print_success "All large files are already being ignored by git."
        exit 0
    fi
    
    print_info "Files to add to .gitignore: ${#files_to_add[@]}"
    
    if [ "$DRY_RUN" = true ]; then
        echo
        print_warning "DRY RUN - No changes will be made to .gitignore"
        print_info "Files that would be added:"
        for file_info in "${files_to_add[@]}"; do
            IFS=':' read -r file_path file_size <<< "$file_info"
            local rel_path
            rel_path=$(realpath --relative-to="$GIT_ROOT" "$file_path" 2>/dev/null || echo "$file_path")
            local size_str
            size_str=$(format_file_size "$file_size")
            echo "  $rel_path ($size_str)"
        done
        exit 0
    fi
    
    # Ask for confirmation
    echo
    read -p "Add ${#files_to_add[@]} files to .gitignore? (y/N): " -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Operation cancelled."
        exit 0
    fi
    
    # Add files to .gitignore
    local added_count=0
    for file_info in "${files_to_add[@]}"; do
        IFS=':' read -r file_path file_size <<< "$file_info"
        if add_to_gitignore "$file_path"; then
            added_count=$((added_count + 1))
            local rel_path
            rel_path=$(realpath --relative-to="$GIT_ROOT" "$file_path" 2>/dev/null || echo "$file_path")
            print_success "Added: $rel_path"
        fi
    done
    
    echo
    print_success "Successfully added $added_count files to .gitignore"
    print_info "Don't forget to commit the updated .gitignore file!"
}

# Run main function with all arguments
main "$@"
