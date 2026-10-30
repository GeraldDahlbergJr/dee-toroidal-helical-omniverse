"""Log an exact, checksummed ZIP fallback when Actions artifacts are unavailable.

Recovery: strip GitHub log timestamps, concatenate lines between the BEGIN/END
markers, base64-decode to ZIP and verify the printed SHA256 before extraction.
Only explicitly selected benchmark outputs are included; no environment secrets.
"""
from pathlib import Path
import base64
import hashlib
import io
import os
import zipfile

HERE = Path(__file__).resolve().parent
PATTERNS = ('*.csv', '*.json', '*.png', '*.log', 'environment.txt')


def bundle():
    files = sorted({p for pattern in PATTERNS for p in HERE.glob(pattern) if p.is_file()})
    if not files:
        raise RuntimeError('No diagnostic files available to preserve')
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, path.name)
    data = buffer.getvalue()
    if len(data) > 8*1024*1024:
        raise RuntimeError('Diagnostic ZIP exceeds 8 MiB log-fallback limit')
    # Check byte-for-byte roundtrip before reporting preservation.
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('ZIP integrity check failed')
        for path in files:
            if archive.read(path.name) != path.read_bytes():
                raise RuntimeError('ZIP byte verification failed')
    return data, [p.name for p in files]


if __name__ == '__main__':
    data, names = bundle()
    digest = hashlib.sha256(data).hexdigest()
    print('Artifact upload unavailable: preserving full diagnostic ZIP in this log.')
    print('DIAGNOSTIC_ZIP_SHA256='+digest)
    print('DIAGNOSTIC_ZIP_BYTES='+str(len(data)))
    print('DIAGNOSTIC_ZIP_FILES='+','.join(names))
    print('::group::Recoverable diagnostic ZIP (base64)')
    print('BEGIN_GEOMETRIC_DIAGNOSTIC_ZIP_BASE64')
    encoded = base64.b64encode(data).decode('ascii')
    for start in range(0, len(encoded), 120):
        print(encoded[start:start+120])
    print('END_GEOMETRIC_DIAGNOSTIC_ZIP_BASE64')
    print('::endgroup::')
    summary = os.environ.get('GITHUB_STEP_SUMMARY')
    if summary:
        with open(summary, 'a') as f:
            f.write('\nArtifact upload failed. Full diagnostics were preserved as a verified '
                    'base64 ZIP in the fallback step log; download the job logs to recover.\n\n')
            f.write(f'ZIP SHA256: `{digest}`; bytes: {len(data)}.\n')
