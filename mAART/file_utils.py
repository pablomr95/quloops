import os
import h5py
import tempfile
import shutil
from filelock import FileLock  # Requires: pip install filelock

def safe_h5_write(filename, mode='w', timeout=10):
    """Safely create/open HDF5 file with locking"""
    lock_path = f"{filename}.lock"
    try:
        with FileLock(lock_path, timeout=timeout):
            # Create parent directories if needed
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            
            # Use temp file then atomic rename
            with tempfile.NamedTemporaryFile(delete=False) as tmp:
                tmp_path = tmp.name
            
            h5f = h5py.File(tmp_path, mode)
            def cleanup():
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
                if os.path.exists(lock_path):
                    os.remove(lock_path)
            
            h5f._tmp_path = tmp_path
            h5f._cleanup = cleanup
            return h5f
    except Exception as e:
        if os.path.exists(lock_path):
            os.remove(lock_path)
        raise RuntimeError(f"Failed to create {filename}: {str(e)}")

def finalize_h5(h5f, target_path):
    """Properly close and move temporary HDF5 file"""
    try:
        h5f.close()
        os.replace(h5f._tmp_path, target_path)
    finally:
        h5f._cleanup()

def validate_hdf5(filename):
    """Thorough HDF5 validation"""
    try:
        with FileLock(f"{filename}.lock", timeout=5):
            with h5py.File(filename, 'r') as f:
                # Check file structure and data integrity
                if not f.keys():
                    return False
                for key in f.keys():
                    if not isinstance(f[key], (h5py.Dataset, h5py.Group)):
                        return False
                return True
    except Exception:
        return False
