import os
import ast
import sys
import importlib
from typing import Set, Dict, List

def get_py_files(directory: str) -> List[str]:
    """Get all .py files excluding temporary/dependency directories"""
    skip_dirs = {"venv", ".venv", "env", "__pycache__", "tests"}
    py_files = []
    for root, dirs, files in os.walk(directory):
        dirs[:] = [d for d in dirs if d not in skip_dirs]
        for file in files:
            if file.endswith(".py"):
                py_files.append(os.path.join(root, file))
    return py_files

def get_local_modules(directory: str) -> Set[str]:
    """Extract all local modules (including submodules)"""
    local_modules = set()
    py_files = get_py_files(directory)

    for file in py_files:
        rel_path = os.path.relpath(file, directory)
        module_path = rel_path.replace(".py", "").replace(os.sep, ".")
        
        if file.endswith("__init__.py"):
            module_path = os.path.dirname(module_path)
        
        parts = module_path.split(".")
        for i in range(len(parts)):
            local_modules.add(".".join(parts[:i+1]))

    return local_modules

def extract_imports(file_path: str, local_modules: Set[str]) -> Set[str]:
    """Extract non-local imports"""
    with open(file_path, "r", encoding="utf-8") as f:
        code = f.read()

    try:
        tree = ast.parse(code)
    except SyntaxError:
        return set()

    imports = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                pkg = alias.name
                if not any(pkg.startswith(local + ".") or pkg == local 
                          for local in local_modules):
                    imports.add(pkg.split(".")[0])

        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                pkg = node.module
                if not any(pkg.startswith(local + ".") or pkg == local 
                          for local in local_modules):
                    imports.add(pkg.split(".")[0])

    return imports

def analyze_dependencies(directory: str) -> Dict[str, Set[str]]:
    """Analyze dependencies across all Python files"""
    py_files = get_py_files(directory)
    local_modules = get_local_modules(directory)
    dependency_map = {}

    for file in py_files:
        imports = extract_imports(file, local_modules)
        dependency_map[file] = imports

    return dependency_map

def print_report(dependency_map: Dict[str, Set[str]]):
    """Generate a dependency report"""
    all_deps = set()
    print("===== Dependency Analysis (Excluding Local Modules) =====")

    for file, deps in dependency_map.items():
        print(f"\n [File] {file}")
        if not deps:
            print("No external dependencies detected")
            continue

        print("Dependencies:")
        for dep in sorted(deps):
            status = "done √" if check_installed(dep) else "Not installed ×"
            print(f"    - {dep} {status}")
            all_deps.add(dep)

    if all_deps:
        print("\n===== Summary =====")
        print("All detected external dependencies:")
        for dep in sorted(all_deps):
            print(f"    - {dep}")

        print("\nCommand to install missing packages:")
        missing_deps = [dep for dep in all_deps if not check_installed(dep)]
        if missing_deps:
            print(f"pip install {' '.join(missing_deps)}")
        else:
            print("All dependencies are installed! √")

def check_installed(package: str) -> bool:
    """Check if a package is installed"""
    try:
        importlib.import_module(package)
        return True
    except ImportError:
        return False

if __name__ == "__main__":
    import argparse

    # 添加当前目录到 sys.path，确保本地模块可被导入
    script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"script")
    print(script_dir)
    sys.path.insert(0, script_dir)

    parser = argparse.ArgumentParser(description="Python dependency analyzer")
    parser.add_argument("directory", help="Project directory", default=".", nargs="?")
    args = parser.parse_args()

    target_dir = os.path.abspath(args.directory)
    if not os.path.isdir(target_dir):
        print("Error: Directory not found")
        exit(1)

    dependency_map = analyze_dependencies(target_dir)
    print_report(dependency_map)
