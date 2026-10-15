"""Independent engine adapter for the admitted flower's existing protocols."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import mujoco
from flexible_bearing import support_protocol


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def nums(values): return ' '.join(format(float(x), '.17g') for x in values)


def load_input(bundle):
    bundle = Path(bundle)
    packet = json.loads((bundle/'table-scene.json').read_text())
    for name, expected in packet['files'].items():
        if sha(bundle/name) != expected: raise ValueError('MATERIAL_IDENTITY_MISMATCH:'+name)
    if packet['classification'] != 'RECONSTRUCTED_METHOD_VALIDATION':
        raise ValueError('INPUT_CLASSIFICATION_MISMATCH')
    return packet


def model_xml(packet, variant, kind, dt):
    if kind not in ('ideal','flexible') or variant not in ('original','optimized'):
        raise ValueError('INVALID_SCENE')
    prop = packet['properties'][variant]
    I = np.asarray(prop['inertia_kg_m2']); c = np.asarray(prop['center_m']); mass = prop['mass_kg']
    if mass <= 0 or I.shape != (3,3) or c.shape != (3,) or not np.isfinite(np.r_[I.ravel(),c,mass]).all() or np.linalg.eigvalsh(I).min() <= 0 or not np.allclose(I,I.T,atol=1e-18):
        raise ValueError('INVALID_MASS_PROPERTIES')
    support = support_protocol(packet['properties']['original']['mass_kg'], packet['properties']['clamp']['mass_kg'])
    root = ET.Element('mujoco', model='admitted-flower-'+kind)
    ET.SubElement(root,'compiler',angle='radian',inertiafromgeom='false')
    opt=ET.SubElement(root,'option',timestep=str(dt),gravity='0 0 -9.81',integrator='RK4')
    ET.SubElement(opt,'flag',energy='enable')
    parent=ET.SubElement(root,'worldbody')
    if kind == 'flexible':
        for j,axis in enumerate(('1 0 0','0 1 0')):
            parent=ET.SubElement(parent,'body',name='carriage'+str(j))
            m=support['carriage_mass_kg']
            ET.SubElement(parent,'inertial',pos='0 0 0',mass=str(m),diaginertia=nums([m*1e-6]*3))
            ET.SubElement(parent,'joint',name='support'+str(j),type='slide',axis=axis,
                stiffness=str(support['stiffness_N_m'][j]),damping=str(support['damping_N_s_m'][j]))
    rotor=ET.SubElement(parent,'body',name='rotor')
    ET.SubElement(rotor,'inertial',pos=nums(c),mass=str(mass),
        fullinertia=nums([I[0,0],I[1,1],I[2,2],I[0,1],I[0,2],I[1,2]]))
    ET.SubElement(rotor,'joint',name='spin',type='hinge',axis='0 0 1',damping='0',frictionloss='0')
    ET.SubElement(rotor,'site',name='axle',pos='0 0 0',size='.001')
    sensors=ET.SubElement(root,'sensor')
    ET.SubElement(sensors,'force',name='bearing_force',site='axle')
    ET.SubElement(sensors,'torque',name='bearing_torque',site='axle')
    return ET.tostring(root,encoding='unicode'),support


def simulate(packet, variant, kind, dt, duration, phase):
    xml,support=model_xml(packet,variant,kind,dt)
    model=mujoco.MjModel.from_xml_string(xml); data=mujoco.MjData(model)
    bid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'rotor')
    basis=np.empty(9); mujoco.mju_quat2Mat(basis,model.body_iquat[bid]); basis=basis.reshape(3,3)
    loaded=basis@np.diag(model.body_inertia[bid])@basis.T
    prop=packet['properties'][variant]
    np.testing.assert_allclose(loaded,prop['inertia_kg_m2'],rtol=1e-12,atol=1e-18)
    np.testing.assert_allclose(model.body_ipos[bid],prop['center_m'],rtol=1e-12,atol=1e-18)
    if model.nu != 0: raise ValueError('UNEXPECTED_DRIVE')
    data.qpos[-1]=phase; data.qvel[-1]=2*math.pi*packet['initial_hz']
    steps=round(duration/dt); records=np.empty((steps+1,12))
    loss=0.; prev_power=0.; damping=np.asarray(support['damping_N_s_m'])
    for i in range(steps+1):
        if i: mujoco.mj_step(model,data)
        mujoco.mj_forward(model,data)
        xy=data.qpos[:2] if kind=='flexible' else np.zeros(2)
        v=data.qvel[:2] if kind=='flexible' else np.zeros(2)
        power=float(np.dot(damping,v*v)) if kind=='flexible' else 0.
        if i: loss += .5*(power+prev_power)*dt
        prev_power=power
        records[i]=[data.time,*xy,data.qpos[-1],data.qvel[-1],*data.sensordata[:3],sum(data.energy),loss,*v]
    if not np.isfinite(records).all(): raise ValueError('NONFINITE_TRAJECTORY')
    energy=records[:,8]-records[0,8]
    # Potential energy may have an arbitrary reference; initial kinetic energy
    # is independently known at zero support displacement.
    c=np.asarray(prop['center_m']); I=np.asarray(prop['inertia_kg_m2']); m=prop['mass_kg']
    e0=.5*(I[2,2]+m*np.dot(c[:2],c[:2]))*(2*math.pi*packet['initial_hz'])**2
    analytic=m*np.linalg.norm(c[:2])*(2*math.pi*packet['initial_hz'])**2
    return dict(records=records,xml=xml,variant=variant,kind=kind,dt=dt,phase=phase,
        loaded_inertia=loaded.tolist(),actuator_count=model.nu,support=support,
        peak_displacement=float(np.linalg.norm(records[:,1:3],axis=1).max()),
        force_mean=float(np.linalg.norm(records[:,5:7],axis=1).mean()),analytic_force=float(analytic),
        energy_drift=float(abs(energy).max()/e0),
        energy_balance_error=float(abs(energy+records[:,9]).max()/e0))


def run(bundle,out):
    packet=load_input(bundle); out.mkdir(parents=True,exist_ok=False)
    scenes=[]
    for kind in ('ideal','flexible'):
        for variant in ('original','optimized'):
            for phase_id in range(8):
                for dt in (.0005,.00025):
                    result=simulate(packet,variant,kind,dt,8,phase_id*math.pi/4)
                    stem=f'{kind}-{variant}-phase{phase_id}-dt{dt}'
                    data=result.pop('records'); xml=result.pop('xml')
                    np.savez_compressed(out/(stem+'.npz'),records=data,
                        columns=['t','support_x','support_y','angle','omega','Fx_local','Fy_local','Fz_local','energy','damping_loss','vx','vy'])
                    (out/(stem+'.xml')).write_text(xml)
                    result.update(file=stem+'.npz',sha256=sha(out/(stem+'.npz')),xml_sha256=sha(out/(stem+'.xml')),phase_index=phase_id)
                    scenes.append(result)
                    print(json.dumps({k:result[k] for k in ('kind','variant','phase_index','dt','peak_displacement','energy_balance_error')}),flush=True)
    # Fail closed on the independent ideal-bearing benchmark. Flexible support
    # accuracy is reported with time-step comparisons, not assumed from pictures.
    for r in scenes:
        if r['kind']=='ideal' and (r['energy_drift']>1e-6 or abs(r['force_mean']-r['analytic_force'])>1e-6):
            raise ValueError('IDEAL_BENCHMARK_FAILED')
    comparisons=[]
    for kind in ('ideal','flexible'):
        for variant in ('original','optimized'):
            for phase in range(8):
                a,b=[next(r for r in scenes if r['kind']==kind and r['variant']==variant and r['phase_index']==phase and r['dt']==dt) for dt in (.0005,.00025)]
                ca=np.load(out/a['file'])['records']; fi=np.load(out/b['file'])['records'][::2]
                comparisons.append(dict(kind=kind,variant=variant,phase=phase,
                    max_position_difference_m=float(np.linalg.norm(ca[:,1:3]-fi[:,1:3],axis=1).max()),
                    max_omega_difference_rad_s=float(abs(ca[:,4]-fi[:,4]).max())))
    summary=dict(schema='flower-mujoco-existing-protocol-v1',mujoco_version=mujoco.__version__,
        input=str(bundle.resolve()),input_sha256=sha(bundle/'table-scene.json'),scenes=scenes,
        timestep_comparisons=comparisons,actual_integration=True,
        scope='IDEAL_BEARING_AND_REGISTERED_LINEAR_SUPPORT_NOT_CALIBRATED_BEARING_OR_HAND')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.bundle,a.output)
