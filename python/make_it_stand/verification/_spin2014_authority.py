"""Independent original-objective equations on canonical rational moments.

Closed-form first derivatives use exact Fraction or explicitly version-locked
GMP rationals, not the solver's AD, Jacobian, multiplier fit or binary64 views.
This private helper is part of the existing
material verifier; numerical qualification never proves geometry, provenance,
optimality, a completed optimization, or a final material seal.
"""
from fractions import Fraction
import math
import re

from .rational_text import fraction_from_rational_text, fraction_to_rational_text


def _read(value):
    return fraction_from_rational_text(value, require_canonical=True)


def _scalar(settings, name, default, *, positive=False, rational_factory=Fraction):
    value = settings.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + ' must be a finite numeric setting')
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(name + ' is outside the supported range')
    return rational_factory(value)


def authority_arithmetic_identity(arithmetic_backend='FRACTION_V1'):
    """Explicit dependency identity; an unavailable requested backend fails."""
    if arithmetic_backend == 'FRACTION_V1':
        return {'backend':'FRACTION_V1','exact':True}
    if arithmetic_backend != 'GMPY2_2_2_1_EXACT_V1':
        raise ValueError('unsupported authority arithmetic backend')
    import gmpy2
    if gmpy2.version() != '2.2.1' or gmpy2.mp_version() != 'GMP 6.3.0':
        raise ValueError('authority GMP dependency version mismatch')
    return {'backend':arithmetic_backend,'gmpy2':gmpy2.version(),
            'gmp':gmpy2.mp_version(),'exact':True}


def _arithmetic(arithmetic_backend):
    authority_arithmetic_identity(arithmetic_backend)
    if arithmetic_backend == 'FRACTION_V1':
        return Fraction, _read, fraction_to_rational_text
    import gmpy2
    def read(value):
        if not isinstance(value,str) or value.count('/') != 1:
            raise ValueError('rational text must contain one slash')
        a,b = value.split('/')
        if not re.fullmatch(r'-?(?:0|[1-9][0-9]*)',a) or \
           not re.fullmatch(r'(?:0|[1-9][0-9]*)',b) or a == '-0':
            raise ValueError('authority value must be a canonical rational')
        numerator,denominator = gmpy2.mpz(a),gmpy2.mpz(b)
        if denominator <= 0:
            raise ValueError('rational denominator must be positive')
        result = gmpy2.mpq(numerator,denominator)
        if result.numerator != numerator or result.denominator != denominator:
            raise ValueError('authority value must be a canonical rational')
        return result
    def write(value):
        return value.numerator.digits(10)+'/'+value.denominator.digits(10)
    return gmpy2.mpq,read,write


def _independent_authority_physics_values(raw_moments, settings, arithmetic):
    """Return original values/first derivatives without any native calls.

    Moment order is M,Lx,Ly,Lz,Qxx,Qyy,Qzz,Qyz,Qxz,Qxy. Units are applied
    before evaluating the formulas; derivatives are with respect to the raw
    input moments. Density is already included and is never applied here.
    """
    rational,read,_ = arithmetic
    names = {'YOYO', 'TOP_INERTIA_ONLY', 'TOP_LOWERING_ONLY',
             'TOP_MASS_AND_LOWERING', 'STATIC_LOWERING'}
    objective = settings.get('objective')
    if objective not in names:
        raise ValueError('unsupported original objective')
    if len(raw_moments) != 10:
        raise ValueError('expected ten authority moments')
    raw = tuple(read(x) for x in raw_moments)
    mass_unit = _scalar(settings,'solver_mass_unit_in_declared_units',1,positive=True,rational_factory=rational)
    length_unit = _scalar(settings,'solver_length_unit_in_declared_units',1,positive=True,rational_factory=rational)
    gamma_i = _scalar(settings,'inertia_weight',1,rational_factory=rational)
    gamma_c = _scalar(settings,'center_weight',1,rational_factory=rational)
    scales = (mass_unit,) + (mass_unit*length_unit,)*3 + (mass_unit*length_unit**2,)*6
    s = tuple(x*d for x,d in zip(raw,scales))
    mass = s[0]
    if mass <= 0:
        raise ValueError('original objective requires positive mass')
    first = s[1:4]
    height = first[2]/mass
    gradient = [rational(0)]*10
    values = [first[0], first[1]]
    rows = [[rational(int(k == i)) for k in range(10)] for i in (1,2)]
    value = rational(0)

    if objective != 'STATIC_LOWERING':
        second_indices = {(0,0):4, (1,1):5, (2,2):6, (1,2):7, (0,2):8, (0,1):9}
        def centered(i,j):
            index = second_indices[tuple(sorted((i,j)))]
            derivative = [rational(0)]*10
            derivative[0] = first[i]*first[j]/mass**2
            derivative[i+1] -= first[j]/mass
            derivative[j+1] -= first[i]/mass
            derivative[index] = 1
            return s[index]-first[i]*first[j]/mass, derivative
        cxx,gxx = centered(0,0)
        cyy,gyy = centered(1,1)
        czz,gzz = centered(2,2)
        cxy,gxy = centered(0,1)
        cxz,gxz = centered(0,2)
        cyz,gyz = centered(1,2)
        ia,ib,ic,ixy = cyy+czz,cxx+czz,cxx+cyy,-cxy
        if abs(ic) <= rational(1e-14):
            raise ValueError('axial central inertia denominator is unqualified')
        ga = [a+b for a,b in zip(gyy,gzz)]
        gb = [a+b for a,b in zip(gxx,gzz)]
        gc = [a+b for a,b in zip(gxx,gyy)]
        gxy = [-x for x in gxy]
        numerator = ia**2+ib**2+2*ixy**2
        value = gamma_i*numerator/ic**2
        gradient = [gamma_i*((2*ia*a+2*ib*b+4*ixy*d)/ic**2-2*numerator*c/ic**3)
                    for a,b,c,d in zip(ga,gb,gc,gxy)]
        if objective == 'YOYO':
            values.append(first[2])
            rows.append([rational(int(k == 3)) for k in range(10)])
        values.extend((-cxz,-cyz))
        rows.extend(([-x for x in gxz],[-x for x in gyz]))
    if objective in ('STATIC_LOWERING','TOP_LOWERING_ONLY'):
        value += gamma_c*height**2
        gradient[0] -= 2*gamma_c*first[2]**2/mass**3
        gradient[3] += 2*gamma_c*first[2]/mass**2
    elif objective == 'TOP_MASS_AND_LOWERING':
        value += gamma_c*first[2]**2
        gradient[3] += 2*gamma_c*first[2]
    return {'objective':value,'objective_gradient':[g*d for g,d in zip(gradient,scales)],
            'constraints':values,'constraint_jacobian':[[g*d for g,d in zip(row,scales)] for row in rows]}


def _serialize_physics(physics,write):
    return {'scope':'INDEPENDENT_NUMERICAL_EQUATIONS_NOT_GEOMETRY_OR_FINAL_SEAL',
            'objective':write(physics['objective']),
            'objective_gradient':[write(g) for g in physics['objective_gradient']],
            'constraints':[write(x) for x in physics['constraints']],
            'constraint_jacobian':[[write(g) for g in row] for row in physics['constraint_jacobian']]}


def independent_authority_physics(raw_moments, settings, *, arithmetic_backend='FRACTION_V1'):
    arithmetic = _arithmetic(arithmetic_backend)
    return _serialize_physics(_independent_authority_physics_values(raw_moments,settings,arithmetic),arithmetic[2])


def check_independent_authority_stationarity(raw_moments, settings, dual_authority,
                                            multipliers, multiplier_mu, *, arithmetic_backend='FRACTION_V1'):
    """Compare exact squared original residual with (10^-12)^2.

    No square root or binary64 roundtrip is needed at the decision boundary.
    Supplied multipliers are checked, not fitted by this verifier. Validating
    their source binding and the geometry is the enclosing material audit's
    responsibility. Passing this equation check alone cannot grant completion.
    """
    rational,read,write = arithmetic = _arithmetic(arithmetic_backend)
    physics = _independent_authority_physics_values(raw_moments,settings,arithmetic)
    if len(dual_authority) != 10 or len(multipliers) != len(physics['constraints']):
        raise ValueError('wrong dual or multiplier dimension')
    dual = tuple(read(x) for x in dual_authority)
    if abs(sum((x*x for x in dual),rational(0))-1) > rational(1,2**1000):
        raise ValueError('authority dual is not unit normalized')
    lagrange = tuple(read(x) for x in multipliers)
    mu = read(multiplier_mu)
    gradient = physics['objective_gradient']
    rows = physics['constraint_jacobian']
    constraints = tuple(physics['constraints'])
    residual = tuple(gradient[k]-sum((row[k]*lam for row,lam in zip(rows,lagrange)),rational(0))+mu*dual[k]
                     for k in range(10)) + constraints
    squared = sum((x*x for x in residual),rational(0))
    within = squared <= rational(1,10**24)
    return {**_serialize_physics(physics,write),
            'raw_equation_components':[write(x) for x in residual],
            'raw_equation_residual_squared':write(squared),
            'raw_equation_within_tolerance':within,
            'positive_multiplier':mu > 0,
            'numerical_stationarity_qualified':within and mu > 0}


def _check_authority_material_moments(point, expected_schema, *, expected_request_identity,
                                    expected_context_identity, arithmetic_backend='FRACTION_V1'):
    """Reintegrate the exported material, without native reconstruction.

    The material region/partition, complete provenance, derivative applicability
    and persistence remain separate gates. These input identities must come from
    the enclosing run; matching strings alone do not authenticate geometry.
    Nothing in this report grants a final seal or a completed optimization.
    """
    authority_arithmetic_identity(arithmetic_backend)
    from .mass_moments import integrate_exact_rational_raw_moments
    for key,expected in (('request_identity',expected_request_identity),
                         ('context_identity',expected_context_identity)):
        if not isinstance(expected,str) or re.fullmatch('[0-9a-f]{64}',expected) is None or point[key] != expected:
            raise ValueError('authority input identity mismatch: '+key)
    if point['schema'] != expected_schema or \
       point['coefficient_policy'] != 'COEFFICIENT_MPFR1024_V1' or point['precision_bits'] != 1024:
        raise ValueError('unsupported authority point policy')
    geometry = tuple(_read(x) for x in point['geometry_authority'])
    dual = tuple(_read(x) for x in point['dual_authority'])
    if len(geometry) != 10 or len(dual) != 10:
        raise ValueError('wrong authority coefficient dimension')
    pairing = tuple(x if i < 7 else 2*x for i,x in enumerate(geometry))
    arithmetic_bound = Fraction(1,2**1000)
    if abs(sum((x*x for x in geometry),Fraction(0))-1) > arithmetic_bound or \
       abs(sum((x*x for x in dual),Fraction(0))-1) > arithmetic_bound or \
       sum((x*y for x,y in zip(pairing,dual)),Fraction(0)) <= 0 or \
       any(abs(pairing[i]*dual[j]-pairing[j]*dual[i]) > arithmetic_bound
           for i in range(10) for j in range(i)):
        raise ValueError('inconsistent authority geometry/dual pair')
    forward = tuple(tuple(_read(x) for x in row) for row in point['model_to_target_authority'])
    inverse = tuple(tuple(_read(x) for x in row) for row in point['target_to_model_authority'])
    if any(len(a) != 3 or any(len(row) != 4 for row in a) for a in (forward,inverse)):
        raise ValueError('invalid authority affine shape')
    for a,b in ((forward,inverse),(inverse,forward)):
        for i in range(3):
            if a[i][3] + sum((a[i][j]*b[j][3] for j in range(3)),Fraction(0)) != 0:
                raise ValueError('affine translations are not inverse')
            for k in range(3):
                if sum((a[i][j]*b[j][k] for j in range(3)),Fraction(0)) != int(i == k):
                    raise ValueError('affine linear maps are not inverse')
    def vertices(field):
        rows = tuple(tuple(_read(x) for x in row) for row in point[field])
        if not rows or any(len(row) != 3 for row in rows):
            raise ValueError('invalid authority vertices: '+field)
        return rows
    def transform(rows,affine):
        return tuple(tuple(affine[i][3]+sum((affine[i][j]*row[j] for j in range(3)),Fraction(0))
                           for i in range(3)) for row in rows)
    source = vertices('source_vertices_target_authority')
    lower = vertices('source_coordinate_lower')
    upper = vertices('source_coordinate_upper')
    if len(source) != len(lower) or len(source) != len(upper):
        raise ValueError('authority source interval count mismatch')
    for coordinates,lo,hi in zip(source,lower,upper):
        for x,a,b in zip(coordinates,lo,hi):
            if not a <= x <= b or b-a > Fraction(1,2**384):
                raise ValueError('authority source coordinate interval invalid')
    outer = transform(vertices('outer_vertices_model_authority'),forward)
    design = transform(vertices('design_vertices_model_authority'),forward)
    if design != vertices('design_vertices_target_authority'):
        raise ValueError('authority design does not match the input affine')
    fixed_density = _scalar(point,'fixed_material_density',None,positive=True)
    design_density = _scalar(point,'design_material_density',None,positive=True)
    om = integrate_exact_rational_raw_moments(outer,point['outer_faces'],arithmetic_backend=arithmetic_backend)
    dm = integrate_exact_rational_raw_moments(design,point['design_faces'],arithmetic_backend=arithmetic_backend)
    fixed = tuple((a-b)*fixed_density for a,b in zip(om,dm))
    if fixed[0] < 0:
        raise ValueError('negative fixed material mass')
    variable = tuple(x*design_density for x in integrate_exact_rational_raw_moments(
        vertices('material_vertices_target_authority'),point['material_faces'],arithmetic_backend=arithmetic_backend))
    total = tuple(a+b for a,b in zip(fixed,variable))
    if variable[0] < 0 or total[0] <= 0:
        raise ValueError('invalid variable or total authority material mass')
    for field,independent in (('fixed_moments_target',fixed),('variable_moments_target',variable),
                              ('total_moments_target',total)):
        if tuple(_read(x) for x in point[field]) != independent:
            raise ValueError('independently integrated authority moments mismatch: '+field)
    write = fraction_to_rational_text
    return dict(scope='INDEPENDENT_MATERIAL_MOMENTS_NOT_GEOMETRY_OR_FINAL_SEAL',
                independent_fixed_moments_target=[write(x) for x in fixed],
                independent_variable_moments_target=[write(x) for x in variable],
                independent_total_moments_target=[write(x) for x in total])


def check_independent_authority_material_mapping(mapping, *, expected_request_identity,
                                               expected_context_identity, arithmetic_backend='FRACTION_V1'):
    """Independent integration of the production mapping; no solver verdict.

    Source authenticity, containment/partition and Jacobian correctness are
    separate gates. This routine deliberately requires no Newton multipliers
    and cannot grant stationarity, optimization or manufacturing approval.
    """
    if mapping.get('scope') != 'MOMENT_MAPPING_NOT_SOLVER_OR_FINAL_SEAL':
        raise ValueError('wrong authority material mapping scope')
    result = _check_authority_material_moments(mapping,'spin2014-authority-material-map-v1',
        expected_request_identity=expected_request_identity,
        expected_context_identity=expected_context_identity,arithmetic_backend=arithmetic_backend)
    if arithmetic_backend!='FRACTION_V1':
        result['arithmetic_identity']=authority_arithmetic_identity(arithmetic_backend)
    return result


def check_independent_authority_point(point, settings, *, expected_request_identity,
                                     expected_context_identity, arithmetic_backend='FRACTION_V1'):
    """Independent material integration and original equations; not a seal."""
    material = _check_authority_material_moments(point,'spin2014-authority-actual-point-v1',
        expected_request_identity=expected_request_identity,
        expected_context_identity=expected_context_identity,arithmetic_backend=arithmetic_backend)
    checked = check_independent_authority_stationarity(material['independent_total_moments_target'],settings,
        point['dual_authority'],point['multipliers'],point['multiplier_mu'],
        arithmetic_backend=arithmetic_backend)
    rational,read,_ = _arithmetic(arithmetic_backend)
    def matches(reported,independent,label):
        if len(reported) != len(independent):
            raise ValueError('authority equation shape mismatch: '+label)
        for a,b in zip(reported,independent):
            exact = read(b)
            if abs(read(a)-exact) > rational(1,10**100)*max(1,abs(exact)):
                raise ValueError('independent authority equation mismatch: '+label)
    matches([point['objective']],[checked['objective']],'objective')
    matches(point['objective_gradient'],checked['objective_gradient'],'gradient')
    matches(point['constraints'],checked['constraints'],'constraints')
    if len(point['constraint_jacobian']) != len(checked['constraint_jacobian']):
        raise ValueError('authority constraint Jacobian shape mismatch')
    for actual,wanted in zip(point['constraint_jacobian'],checked['constraint_jacobian']):
        matches(actual,wanted,'constraint Jacobian')
    raw = read(point['raw_equation_residual'])
    if raw < 0:
        raise ValueError('authority raw residual cannot be negative')
    squared = read(checked['raw_equation_residual_squared'])
    if abs(raw**2-squared) > rational(1,10**100)*max(1,abs(squared)):
        raise ValueError('independent original-equation residual mismatch')
    return {**checked, **material,
            'scope':'INDEPENDENT_MATERIAL_MOMENTS_AND_EQUATIONS_NOT_GEOMETRY_OR_FINAL_SEAL'}


def check_authority_initial_objective(attempt, settings):
    """Check the original objective on recorded density-weighted initial moments.

    The enclosing run must independently bind/reconstruct this initial material;
    this arithmetic check alone is not an initial geometry/provenance proof.
    """
    moments = attempt.get('initial_total_moments')
    identity = attempt.get('initial_material_identity')
    if not isinstance(moments,(tuple,list)) or not isinstance(identity,str) or \
       re.fullmatch('[0-9a-f]{64}',identity) is None:
        raise ValueError('initial authority material evidence missing')
    physics = independent_authority_physics(moments,settings)
    actual = _read(attempt['initial_objective'])
    independent = _read(physics['objective'])
    if abs(actual-independent) > Fraction(1,10**100)*max(1,abs(independent)):
        raise ValueError('initial original objective disagrees with recorded material moments')
    return dict(scope='INITIAL_OBJECTIVE_ARITHMETIC_NOT_GEOMETRY_PROOF',
                independent_initial_objective=physics['objective'])


def check_authority_event_journal(path, *, run, expected_request_identity,
                                  expected_context_identity, expected_objective_identity):
    """Read the durable stream and compare with a completed native return.

    A truncated prefix remains useful evidence but cannot pass this completed
    stream check. Equality is observation evidence, never material admission.
    Expected identities come from the enclosing frozen-input audit.
    """
    import hashlib
    import json
    expected = dict(request_identity=expected_request_identity,
                    context_identity=expected_context_identity,
                    objective_identity=expected_objective_identity)
    for key,value in expected.items():
        if not isinstance(value,str) or re.fullmatch('[0-9a-f]{64}',value) is None:
            raise ValueError('invalid expected authority event identity')
        if key != 'objective_identity' and run.get(key) != value:
            raise ValueError('authority run/event input identity mismatch')
    events = run.get('events')
    if not isinstance(events,list) or not events or events[0].get('type') != 'attempt_start' or \
       events[-1].get('type') != 'run_terminal':
        raise ValueError('authority returned stream is not complete')
    def unique_object(pairs):
        result = {}
        for key,value in pairs:
            if key in result:
                raise ValueError('duplicate authority journal JSON key')
            result[key] = value
        return result
    def invalid_constant(value):
        raise ValueError('nonfinite authority journal JSON value')
    count = 0
    digest = hashlib.sha256()
    with open(path,'r',encoding='utf-8') as stream:
        for line in stream:
            digest.update(line.encode('utf-8'))
            if not line.endswith('\n') or count >= len(events):
                raise ValueError('truncated or extra authority journal record')
            row = json.loads(line,object_pairs_hook=unique_object,parse_constant=invalid_constant)
            if not isinstance(row,dict) or row.get('schema') != 'spin2014-solver-event-journal-v1' or \
               type(row.get('sequence')) is not int or row['sequence'] != count:
                raise ValueError('authority journal sequence/schema mismatch')
            event = row.get('event')
            if not isinstance(event,dict) or type(event.get('sequence')) is not int or \
               event['sequence'] != count or event != events[count]:
                raise ValueError('authority journal differs from returned native stream')
            if any(event.get(key) != value for key,value in expected.items()):
                raise ValueError('authority journal effective input identity mismatch')
            count += 1
    if count != len(events):
        raise ValueError('missing authority journal records')
    return dict(scope='DURABLE_OBSERVATION_NOT_FINAL_ADMISSION',
                completed_stream_matches=True,records=count,
                journal_sha256=digest.hexdigest(),sealed_artifact=False)


def check_authority_restart_ledger(run, *, objective_settings, expected_initial_sequence, evaluation_budget_per_attempt,
                                   expected_request_identity, expected_context_identity):
    """Check actual attempt accounting; NOT a numerical or provenance proof.

    Expected starts must be derived from the frozen request/seed, not this run.
    Full candidate eligibility still needs independent point/material checks,
    actual material change, source reconstruction and complete final run admission.
    """
    def integer(value,label,maximum):
        if type(value) is not int or not 0 <= value <= maximum:
            raise ValueError('invalid authority ledger integer: '+label)
        return value
    def identity(value,label):
        if not isinstance(value,str) or re.fullmatch('[0-9a-f]{64}',value) is None:
            raise ValueError('invalid authority ledger identity: '+label)
    budget = integer(evaluation_budget_per_attempt,'budget',60)
    if budget == 0 or not 1 <= len(expected_initial_sequence) <= 8:
        raise ValueError('invalid authority attempt protocol')
    if run['schema'] != 'spin2014-authority-restarts-v1':
        raise ValueError('wrong authority restart schema')
    for key,expected in (('request_identity',expected_request_identity),
                         ('context_identity',expected_context_identity)):
        identity(expected,key)
        if run[key] != expected:
            raise ValueError('authority ledger input mismatch: '+key)
    scope = ('CONFIGURED_FULL_ATTEMPT_BUDGET' if budget == 60
             else 'DIAGNOSTIC_REDUCED_ATTEMPT_BUDGET')
    if run['protocol_scope'] != scope:
        raise ValueError('authority ledger protocol scope mismatch')
    attempts = run['attempts']
    index = run['stationary_attempt_index']
    if not isinstance(attempts,list) or not 1 <= len(attempts) <= len(expected_initial_sequence):
        raise ValueError('invalid authority attempt count')
    if index is None:
        if run['stationary_point'] is not None or run['status'] != 'restart_budget_exhausted' or \
           len(attempts) != len(expected_initial_sequence):
            raise ValueError('failed authority attempts missing or candidate status inconsistent')
    elif type(index) is not int or index != len(attempts)-1 or run['status'] != 'stationary_candidate':
        raise ValueError('authority stationary attempt index inconsistent')
    total = 0
    for number,attempt in enumerate(attempts):
        expected = expected_initial_sequence[number]
        for a,b in (('initial_geometry_authority','geometry_authority'),
                    ('initial_dual_authority','dual_authority')):
            if attempt[a] != expected[b] or len(attempt[a]) != 10:
                raise ValueError('authority initial sequence mismatch')
            for value in attempt[a]:
                _read(value)
        count = integer(attempt['function_evaluations'],'evaluations',budget)
        if attempt['initial_state_identity']:
            identity(attempt['initial_state_identity'],'initial state')
            check_authority_initial_objective(attempt,objective_settings)
        elif attempt.get('initial_total_moments') is not None or attempt.get('initial_objective') or \
             attempt.get('initial_material_identity'):
            raise ValueError('failed initial construction carries unbound material evidence')
        trials = attempt['trials']
        if count != 1+len(trials):
            raise ValueError('authority initial/trial evaluation accounting mismatch')
        accepted = []
        for trial in trials:
            if type(trial['accepted']) is not bool:
                raise ValueError('authority accepted flag must be Boolean')
            arc = _read(trial['arc_length'])
            if not 0 < arc <= Fraction(1,2):
                raise ValueError('authority trial arc outside frozen protocol')
            if trial['raw_equation_residual'] and _read(trial['raw_equation_residual']) < 0:
                raise ValueError('negative authority trial residual')
            if trial['accepted']:
                if trial['reason'] not in ('residual_decreased','topology_transition_relinearized') or \
                   not trial['raw_equation_residual']:
                    raise ValueError('authority accepted trial lacks actual decrease evidence')
                accepted.append(trial)
        updates = integer(attempt['accepted_updates'],'accepted updates',count)
        if updates != len(accepted) or updates != len(attempt['accepted_state_identities']) or \
           updates > len(attempt['step_identities']):
            raise ValueError('authority accepted update accounting mismatch')
        for value in attempt['step_identities']+attempt['accepted_state_identities']:
            identity(value,'step/state')
        if updates and attempt['final_state_identity'] != attempt['accepted_state_identities'][-1]:
            raise ValueError('authority final state does not match the last accepted update')
        candidate = attempt['status'] == 'stationary_candidate'
        if candidate != (index == number):
            raise ValueError('authority attempt candidate label mismatch')
        if candidate:
            point = run['stationary_point']
            if not isinstance(point,dict) or point['state_scientific_hash'] != attempt['final_state_identity'] or \
               point['request_identity'] != expected_request_identity or \
               point['context_identity'] != expected_context_identity or \
               point['raw_equation_residual'] != attempt['final_raw_residual'] or \
               not 0 <= _read(attempt['final_raw_residual']) <= Fraction(1,10**12) or \
               _read(point['multiplier_mu']) <= 0:
                raise ValueError('authority terminal point does not match the candidate attempt')
        total += count
    if type(run['function_evaluations']) is not int or run['function_evaluations'] != total:
        raise ValueError('authority restart total evaluation accounting mismatch')
    return dict(scope='STRUCTURAL_RESTART_LEDGER_NOT_FINAL_ADMISSION',function_evaluations=total,
                full_run_candidate=budget==60 and index is not None and attempts[index]['accepted_updates']>0,
                sealed_artifact=False)


def require_authority_initial_screen(point):
    """Exact screen on bound production values; derivative proof is separate."""
    if _read(point['multiplier_mu'])<=0:
        raise ValueError('nonpositive initial multiplier')
    matrix=point['native_moment_jacobian']
    if not isinstance(matrix,(list,tuple)) or len(matrix)!=10 or any(
            not isinstance(row,(list,tuple)) or len(row)!=10 for row in matrix):
        raise ValueError('invalid initial moment Jacobian dimension')
    values=tuple(_read(value) for row in matrix for value in row)
    if not any(value!=0 for value in values):
        raise ValueError('zero initial moment sensitivity')


def check_authority_complete_run_evidence(run, *, settings, expected_initial_sequence,
        journal_path, expected_request_identity, expected_context_identity,
        expected_objective_identity, initial_point, arithmetic_backend='GMPY2_2_2_1_EXACT_V1'):
    """Independent complete-run numerical gate, NOT geometry or persistence seal.

    Expected starts/identities must be reconstructed from frozen input by the
    enclosing caller. Initial point must be actual reconstructed material,
    never a hand-filled moment vector or a diagnostic replacing a restart run.
    """
    if run.get('schema')!='spin2014-authority-restarts-v1' or run.get('protocol_scope')!='CONFIGURED_FULL_ATTEMPT_BUDGET':
        raise ValueError('complete full-budget run required')
    if not isinstance(initial_point,dict) or initial_point.get('schema')!='spin2014-authority-actual-point-v1':
        raise ValueError('initial actual material required')
    ledger=check_authority_restart_ledger(run,objective_settings=settings,
        expected_initial_sequence=expected_initial_sequence,evaluation_budget_per_attempt=60,
        expected_request_identity=expected_request_identity,expected_context_identity=expected_context_identity)
    if not ledger['full_run_candidate']:
        raise ValueError('actual accepted optimization updates required')
    journal=check_authority_event_journal(journal_path,run=run,
        expected_request_identity=expected_request_identity,expected_context_identity=expected_context_identity,
        expected_objective_identity=expected_objective_identity)
    attempt=run['attempts'][run['stationary_attempt_index']]
    point=run['stationary_point']
    for key,value in (('state_scientific_hash',attempt['initial_state_identity']),
                      ('geometry_authority',attempt['initial_geometry_authority']),
                      ('dual_authority',attempt['initial_dual_authority']),
                      ('material_scientific_hash',attempt['initial_material_identity']),
                      ('multiplier_mu',attempt['initial_mu']),
                      ('raw_equation_residual',attempt['initial_raw_residual']),
                      ('total_moments_target',attempt['initial_total_moments'])):
        if initial_point.get(key)!=value:
            raise ValueError('initial actual material binding mismatch: '+key)
    for field in ('outer_vertices_model_authority','outer_faces','design_vertices_model_authority',
                  'design_faces','model_to_target_authority','target_to_model_authority',
                  'fixed_material_density','design_material_density','fixed_moments_target'):
        if initial_point.get(field)!=point.get(field):
            raise ValueError('fixed input changed during run: '+field)
    if initial_point.get('objective_identity')!=expected_objective_identity or point.get('objective_identity')!=expected_objective_identity:
        raise ValueError('objective identity mismatch')
    initial=check_independent_authority_point(initial_point,settings,
        expected_request_identity=expected_request_identity,expected_context_identity=expected_context_identity,
        arithmetic_backend=arithmetic_backend)
    require_authority_initial_screen(initial_point)
    terminal=check_independent_authority_point(point,settings,
        expected_request_identity=expected_request_identity,expected_context_identity=expected_context_identity,
        arithmetic_backend=arithmetic_backend)
    if not terminal['numerical_stationarity_qualified']:
        raise ValueError('independent terminal stationarity failed')
    change=check_authority_material_change(attempt,point)
    if not change['material_change_witness']:
        raise ValueError('actual material change not witnessed')
    return dict(scope='COMPLETE_RUN_NUMERICAL_GATE_NOT_GEOMETRY_OR_PERSISTENCE_SEAL',
        ledger=ledger,journal=journal,initial=initial,terminal=terminal,material_change=change,
        eligible_for_material_preparation=True,sealed_artifact=False,manufacturing_approved=False)


def check_authority_material_change(attempt, point):
    """Sufficient moment witness, conditional on immutable geometry/density binding.

    Hash changes alone may reflect coefficients or derivatives rather than
    material. Equal moments do not prove equal materials; they yield UNKNOWN.
    The enclosing admission must authenticate initial/final geometry and density.
    """
    updates=attempt['accepted_updates']
    if type(updates) is not int or updates<0:
        raise ValueError('invalid accepted update count')
    before=attempt['initial_total_moments']; after=point['total_moments_target']
    if not isinstance(before,(list,tuple)) or not isinstance(after,(list,tuple)) or len(before)!=10 or len(after)!=10:
        raise ValueError('material change requires ten initial and final moments')
    difference=tuple(_read(b)-_read(a) for a,b in zip(before,after))
    witness=updates>0 and any(difference)
    return dict(scope='MOMENT_CHANGE_CONDITIONAL_ON_FIXED_INPUT_BINDING_NOT_FINAL_ADMISSION',
        material_change_witness=bool(witness),
        status='CHANGED_MOMENT_WITNESS' if witness else 'UNKNOWN',
        moment_difference=[fraction_to_rational_text(x) for x in difference],sealed_artifact=False)


def prepare_authority_material_boundaries(point, settings, *, expected_request_identity,
                                          expected_context_identity, arithmetic_backend='FRACTION_V1'):
    """Prepare the existing six-boundary representation from the OWNED point.

    No coefficient display, source triangulation, optimization or opposite-q
    evaluation is used. Exact DESIGN-minus-saved-retained is the sole cavity
    construction. This can inspect a nonstationary diagnostic point; its
    numerical result is never promoted to completed optimization or sealing.
    """
    authority_arithmetic_identity(arithmetic_backend)
    from .spin_common import native, require, digest
    from .spin_problem import MeshSnapshot
    from .faster_material import _nested_material_boundary, _canonical_exact_boolean_boundary
    from .mass_moments import integrate_exact_rational_raw_moments
    numerical = check_independent_authority_point(point,settings,
        expected_request_identity=expected_request_identity,
        expected_context_identity=expected_context_identity,arithmetic_backend=arithmetic_backend)
    outer = MeshSnapshot.from_rationals(point['outer_vertices_model_authority'],point['outer_faces'])
    design = MeshSnapshot.from_rationals(point['design_vertices_model_authority'],point['design_faces'])
    inverse = tuple(tuple(_read(x) for x in row) for row in point['target_to_model_authority'])
    model = []
    for row in point['material_vertices_target_authority']:
        source = tuple(_read(x) for x in row)
        model.append([fraction_to_rational_text(inverse[i][3]+sum(
            (inverse[i][j]*source[j] for j in range(3)),Fraction(0))) for i in range(3)])
    retained = _canonical_exact_boolean_boundary(model,point['material_faces'])
    n = native()
    for one,two,label in ((design,outer,'design'),(retained,design,'retained')):
        check = n.check_exact_material_subset(one.exact_region(),two.exact_region())
        require(check['completed'] and check['a_subset_b'],label+' authority containment failed')
    exact = n.libigl_exact_boolean_construction_from_rationals(
        design.rational_vertices,design.arrays()[1],
        retained.rational_vertices,retained.arrays()[1],'difference')
    require(exact['succeeded'],'authority same-retained complement construction failed')
    cavity = (_canonical_exact_boolean_boundary(exact['vertices'],exact['faces'])
              if len(exact['faces']) else MeshSnapshot(None))
    fixed = _nested_material_boundary(outer,design)
    material = _nested_material_boundary(outer,cavity)
    regions = dict(outer=outer,fixed=fixed,design=design,retained=retained,cavity=cavity,material=material)
    for name,mesh in regions.items():
        # A negative cavity shell is reversed ONLY on the checker's temporary
        # positive component copy before subtraction. It is not an input repair.
        # Reject actual input winding errors using the original sign/depth proof.
        facts = mesh.exact_region().facts
        require(len(facts['original_signs']) == len(facts['depths']) and all(
            sign == (1 if depth % 2 == 0 else -1)
            for sign,depth in zip(facts['original_signs'],facts['depths'])),
            name+' authority original winding violates nesting')
    def partition(whole,parts,label):
        facts = n.check_exact_material_partition(whole.exact_region(),[part.exact_region() for part in parts])
        require(facts['completed'] and facts['complete_partition'] and
                facts['no_positive_volume_overlap'] and all(facts['part_subset_whole']),
                label+' authority partition failed')
        return facts
    partitions = dict(design_partition=partition(design,[retained,cavity],'design'),
                      whole_partition=partition(outer,[fixed,retained,cavity],'whole'),
                      material_partition=partition(material,[fixed,retained],'material'))
    def moments(mesh,density):
        if mesh.explicit_empty:
            return (Fraction(0),)*10
        return tuple(x*density for x in integrate_exact_rational_raw_moments(
            mesh.rational_vertices,mesh.arrays()[1].tolist(),arithmetic_backend=arithmetic_backend))
    fixed_model = moments(fixed,_scalar(point,'fixed_material_density',None,positive=True))
    retained_model = moments(retained,_scalar(point,'design_material_density',None,positive=True))
    total_model = tuple(a+b for a,b in zip(fixed_model,retained_model))
    write = fraction_to_rational_text
    audit = dict(numerical,**partitions,
        schema='spin2014-authority-material-boundaries-v1',
        scope='EXACT_REGION_PARTITION_AND_NUMERICAL_POINT_NOT_FINAL_SEAL',
        source_point_digest=digest(point),
        complement_construction='DESIGN_MINUS_SAVED_RETAINED_EXACT_V1',
        fixed_moments_model=[write(x) for x in fixed_model],
        retained_moments_model=[write(x) for x in retained_model],
        total_moments_model=[write(x) for x in total_model],
        manufacturing_approved=False,sealed_artifact=False)
    return regions,audit
