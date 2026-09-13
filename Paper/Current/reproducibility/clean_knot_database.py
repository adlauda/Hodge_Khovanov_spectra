"""Create a unique-bidegree working copy; never modify the published snapshot.

Keep the highest-ID (last-inserted) row for each (knot_id,h,q), but only after
checking its decoded Betti number against the stored completed polynomial.
This reconciles stored records; it is not an independent homology calculation.
"""
import argparse
import hashlib
import json
import re
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT.parents[1] / 'Backup/2026-09-08-R4/reproducibility/snapshots'
SOURCE_COMMIT = 'ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea'
SOURCE_SHA256 = '918832fc674a8b92aeef88f464bad18c4a21478ef90fefa2e57f52aed1d1e21a'

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def decode(value):
    if isinstance(value, bytes):
        if len(value) != 8:
            raise ValueError('Unexpected Betti blob length')
        return int.from_bytes(value, 'little', signed=True)
    return value

def polynomial(text):
    if text is None:
        raise ValueError('Incomplete knot record')
    result = {}
    for term in text.replace(' ', '').split('+'):
        if term in ('', '0'):
            continue
        coeff = re.match(r'^\d*', term).group()
        body = term[len(coeff):]
        tokens = re.findall(r'([qt])(?:\^(-?\d+))?', body)
        if re.sub(r'[qt](?:\^-?\d+)?', '', body) or len({x[0] for x in tokens}) != len(tokens):
            raise ValueError(f'Unrecognized polynomial term: {term}')
        powers = {letter:int(exponent or 1) for letter,exponent in tokens}
        key = (powers.get('t',0), powers.get('q',0))
        if key in result:
            raise ValueError(f'Repeated polynomial monomial: {key}')
        result[key] = int(coeff or 1)
    return result

def check_canonical_betti(db):
    expected = {row[0]:polynomial(row[1]) for row in db.execute('SELECT id,khovanov_polynomial FROM knots')}
    actual = {key:{} for key in expected}
    mismatches = []
    for knot_id,h,q,betti in db.execute('''SELECT knot_id,h,q,betti FROM bidegrees
            WHERE id IN (SELECT MAX(id) FROM bidegrees GROUP BY knot_id,h,q)'''):
        value = decode(betti)
        if value != expected[knot_id].get((h,q),0):
            mismatches.append((knot_id,h,q,value,expected[knot_id].get((h,q),0)))
        if value:
            actual[knot_id][(h,q)] = value
    if mismatches or actual != expected:
        raise ValueError(f'Last-inserted rows disagree with completed polynomials: {mismatches[:10]}')

def main(args):
    source = args.source
    if source is None:
        source = next((path for path in (ROOT/'snapshots/knot_research.db',
                                        ARCHIVE/'knot_research.db') if path.is_file()), None)
        if source is None:
            raise ValueError(f'Pass --source with the original knot_research.db from commit {SOURCE_COMMIT}. '
                             'The current repository database is already cleaned; it is not a cleanup source.')
    source, output, report_path = source.resolve(), args.output.resolve(), args.report.resolve()
    if source == output or output.exists() or report_path.exists():
        raise ValueError('Source and output must differ; output and report must not already exist')
    if digest(source) != SOURCE_SHA256:
        raise ValueError('Unexpected published snapshot checksum')
    report = {'source_commit':SOURCE_COMMIT, 'source_sha256':SOURCE_SHA256,
              'created_utc':datetime.now(timezone.utc).isoformat(),
              'selection_rule':'Keep the highest-ID row per knot and bidegree; verify all retained Betti numbers against the stored completed polynomial.',
              'mathematical_status':'Storage consistency only; no independent spectral or homology certification.',
              'source_unchanged':False, 'conflicting_betti_groups':[], 'removed_row_to_retained_row':[]}
    with sqlite3.connect(source.as_uri()+'?mode=ro',uri=True) as src:
        src.row_factory = sqlite3.Row
        if src.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
            raise ValueError('Source integrity check failed')
        check_canonical_betti(src)
        report['knots'] = src.execute('SELECT COUNT(*) FROM knots').fetchone()[0]
        report['source_bidegree_rows'] = src.execute('SELECT COUNT(*) FROM bidegrees').fetchone()[0]
        repeated = src.execute('''SELECT knot_id,h,q,COUNT(*) n,MAX(id) keep FROM bidegrees
            GROUP BY knot_id,h,q HAVING COUNT(*)>1''').fetchall()
        report['duplicate_groups'] = len(repeated)
        report['knots_with_duplicates'] = len({row['knot_id'] for row in repeated})
        report['multiplicities'] = dict(Counter(row['n'] for row in repeated))
        names = dict(src.execute('SELECT id,name FROM knots'))
        max_spread = 0.0
        for group in repeated:
            rows = src.execute('SELECT * FROM bidegrees WHERE knot_id=? AND h=? AND q=? ORDER BY id',
                               (group['knot_id'],group['h'],group['q'])).fetchall()
            if len({row['dimension'] for row in rows}) != 1:
                raise ValueError('Conflicting chain dimensions; cleanup requires review')
            gaps = [row['smallest_nonzero'] for row in rows]
            if any(g is None for g in gaps) and not all(g is None for g in gaps):
                raise ValueError('Conflicting null and nonnull gap records')
            if gaps[0] is not None:
                spread = max(gaps)-min(gaps)
                max_spread = max(max_spread,spread)
                if min(gaps) <= 0 or spread > 1e-8:
                    raise ValueError('Spectral disagreement requires review')
            report['removed_row_to_retained_row'].extend([row['id'],group['keep']] for row in rows[:-1])
            if len({decode(row['betti']) for row in rows}) > 1:
                report['conflicting_betti_groups'].append({
                    'name':names[group['knot_id']], 'h':group['h'], 'q':group['q'],
                    'retained_row_id':group['keep'],
                    'rows':[{key:decode(row[key]) for key in row.keys()} for row in rows]})
        report['max_duplicate_gap_spread'] = max_spread
        report['conflicting_betti_knots'] = len({r['name'] for r in report['conflicting_betti_groups']})
        output.parent.mkdir(parents=True,exist_ok=True)
        with sqlite3.connect(output.as_uri(),uri=True) as dst:
            src.backup(dst)
            dst.execute('BEGIN IMMEDIATE')
            dst.execute('DELETE FROM bidegrees WHERE id NOT IN (SELECT MAX(id) FROM bidegrees GROUP BY knot_id,h,q)')
            dst.execute('CREATE UNIQUE INDEX unique_knot_bidegree ON bidegrees(knot_id,h,q)')
            # The composite index also supports lookups on its knot_id prefix.
            dst.execute('DROP INDEX idx_bidegrees_knot_id')
            dst.commit()
            check_canonical_betti(dst)
            report['retained_bidegree_rows'] = dst.execute('SELECT COUNT(*) FROM bidegrees').fetchone()[0]
            report['remaining_duplicate_groups'] = dst.execute('''SELECT COUNT(*) FROM (
                SELECT knot_id,h,q FROM bidegrees GROUP BY knot_id,h,q HAVING COUNT(*)>1)''').fetchone()[0]
            report['integrity_check'] = dst.execute('PRAGMA integrity_check').fetchone()[0]
            report['foreign_key_violations'] = dst.execute('PRAGMA foreign_key_check').fetchall()
            dst.execute('ATTACH DATABASE ? AS original',(source.as_uri()+'?mode=ro',))
            report['retained_rows_differing_from_source'] = dst.execute('SELECT COUNT(*) FROM (SELECT * FROM bidegrees EXCEPT SELECT * FROM original.bidegrees)').fetchone()[0]
            for table in ('knots','merged_knot_results'):
                for left,right in ((table,'original.'+table),('original.'+table,table)):
                    if dst.execute(f'SELECT COUNT(*) FROM (SELECT * FROM {left} EXCEPT SELECT * FROM {right})').fetchone()[0]:
                        raise ValueError(f'Unexpected change to {table}')
            report['max_per_diagram_gap_change'] = dst.execute('''SELECT MAX(ABS(a.gap-b.gap)) FROM
                (SELECT knot_id,MIN(smallest_nonzero) gap FROM bidegrees GROUP BY knot_id) a JOIN
                (SELECT knot_id,MIN(smallest_nonzero) gap FROM original.bidegrees GROUP BY knot_id) b USING(knot_id)''').fetchone()[0]
            dst.execute('DETACH DATABASE original')
            dst.execute('VACUUM')
        report['removed_rows'] = report['source_bidegree_rows']-report['retained_bidegree_rows']
    report['output_sha256'] = digest(output)
    report['source_unchanged'] = digest(source) == SOURCE_SHA256
    assert report['source_unchanged'] and report['integrity_check'] == 'ok'
    assert report['remaining_duplicate_groups'] == 0 and report['retained_rows_differing_from_source'] == 0
    assert not report['foreign_key_violations']
    report_path.parent.mkdir(parents=True,exist_ok=True)
    report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('removed_row_to_retained_row','conflicting_betti_groups')},indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,
                        help='Original database from SOURCE_COMMIT; defaults to a local archived snapshot, never the cleaned repository database.')
    parser.add_argument('--output',type=Path,default=ROOT/'data/knot_research.db')
    parser.add_argument('--report',type=Path,default=ROOT/'DATABASE-CLEANUP.json')
    main(parser.parse_args())
