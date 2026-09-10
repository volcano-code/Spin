"""Scoped actual X11 observation; no virtual desktop or scientific file writes.

Use system Python (its Pillow supports XCB). Window ID is independently observed
with xwininfo; title must match this task before capture or pointer input.
"""
import argparse
import ctypes as C
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from PIL import ImageGrab

p = argparse.ArgumentParser()
p.add_argument('action', choices=['capture', 'click', 'snapshot', 'activate'])
p.add_argument('--window')
p.add_argument('--name', required=True)
p.add_argument('--x', type=int)
p.add_argument('--y', type=int)
a = p.parse_args()
assert re.fullmatch(r'[a-z0-9-]+', a.name)
out = Path('build/task-7v1-desktop')
out.mkdir(exist_ok=True)
if a.action == 'snapshot':
    root = Path('build/task-7v1-source-run-1/sealed')
    data = {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in sorted(root.rglob('*')) if f.is_file()}
    target = out/(a.name+'.json')
    with target.open('x') as stream:
        json.dump(data, stream, sort_keys=True, indent=2)
    print(target, len(data), 'sealed files')
else:
    assert a.window and re.fullmatch(r'0x[0-9a-fA-F]+', a.window)
    info = subprocess.check_output(['xwininfo','-display',':1','-id',a.window], text=True)
    assert '"Sealed spin | HOLLOW_ONLY | PARTIAL_BLOCKED"' in info
    assert 'Map State: IsViewable' in info
    def number(label):
        return int(re.search(re.escape(label)+r':\s+(-?\d+)', info).group(1))
    x,y,w,h = [number(s) for s in ('Absolute upper-left X','Absolute upper-left Y','Width','Height')]
    if a.action == 'activate':
        class Data(C.Union):
            _fields_ = [('b',C.c_char*20),('s',C.c_short*10),('l',C.c_long*5)]
        class Message(C.Structure):
            _fields_ = [('type',C.c_int),('serial',C.c_ulong),('send_event',C.c_int),
                        ('display',C.c_void_p),('window',C.c_ulong),('message_type',C.c_ulong),
                        ('format',C.c_int),('data',Data)]
        class Event(C.Union):
            _fields_ = [('message',Message),('pad',C.c_long*24)]
        X=C.CDLL('libX11.so.6')
        X.XOpenDisplay.argtypes=[C.c_char_p]; X.XOpenDisplay.restype=C.c_void_p
        X.XDefaultRootWindow.argtypes=[C.c_void_p]; X.XDefaultRootWindow.restype=C.c_ulong
        X.XInternAtom.argtypes=[C.c_void_p,C.c_char_p,C.c_int]; X.XInternAtom.restype=C.c_ulong
        X.XSendEvent.argtypes=[C.c_void_p,C.c_ulong,C.c_int,C.c_long,C.POINTER(Event)]
        X.XFlush.argtypes=[C.c_void_p]; X.XCloseDisplay.argtypes=[C.c_void_p]
        d=X.XOpenDisplay(b':1'); assert d
        try:
            event=Event(); event.message.type=33; event.message.display=d
            event.message.window=int(a.window,16); event.message.format=32
            event.message.message_type=X.XInternAtom(d,b'_NET_ACTIVE_WINDOW',False)
            event.message.data.l[0]=2
            assert X.XSendEvent(d,X.XDefaultRootWindow(d),False,(1<<20)|(1<<19),C.byref(event))
            X.XFlush(d)
        finally:
            X.XCloseDisplay(d)
        print('requested activation only of named viewer',a.window)
    elif a.action == 'capture':
        active = subprocess.check_output(['xprop','-display',':1','-root','_NET_ACTIVE_WINDOW'],text=True)
        assert int(re.search(r'0x[0-9a-fA-F]+',active).group(),16) == int(a.window,16), active
        target = out/(a.name+'.png')
        assert not target.exists()
        ImageGrab.grab(bbox=(x,y,x+w,y+h), xdisplay=':1').save(target)
        with (out/(a.name+'.json')).open('x') as stream:
            json.dump({'window':a.window,'active_window':active,'xwininfo':info,'image_sha256':hashlib.sha256(target.read_bytes()).hexdigest()},stream,indent=2)
        print(target, (x,y,w,h))
    else:
        assert a.x is not None and a.y is not None and 0 <= a.x < w and 0 <= a.y < h
        active = subprocess.check_output(['xprop','-display',':1','-root','_NET_ACTIVE_WINDOW'],text=True)
        assert int(re.search(r'0x[0-9a-fA-F]+',active).group(),16) == int(a.window,16), active
        X=C.CDLL('libX11.so.6'); T=C.CDLL('libXtst.so.6')
        X.XOpenDisplay.argtypes=[C.c_char_p]; X.XOpenDisplay.restype=C.c_void_p
        X.XRaiseWindow.argtypes=[C.c_void_p,C.c_ulong]
        X.XFlush.argtypes=[C.c_void_p]; X.XCloseDisplay.argtypes=[C.c_void_p]
        T.XTestFakeMotionEvent.argtypes=[C.c_void_p,C.c_int,C.c_int,C.c_int,C.c_ulong]
        T.XTestFakeButtonEvent.argtypes=[C.c_void_p,C.c_uint,C.c_int,C.c_ulong]
        d=X.XOpenDisplay(b':1'); assert d
        try:
            X.XRaiseWindow(d,int(a.window,16)); X.XFlush(d)
            T.XTestFakeMotionEvent(d,-1,x+a.x,y+a.y,0); X.XFlush(d)
            time.sleep(.15)
            T.XTestFakeButtonEvent(d,1,1,0); X.XFlush(d)
            time.sleep(.12)
            T.XTestFakeButtonEvent(d,1,0,0); X.XFlush(d)
        finally:
            X.XCloseDisplay(d)
        with (out/(a.name+'.json')).open('x') as stream:
            json.dump({'window':a.window,'action':'click','relative_xy':[a.x,a.y],'xwininfo':info},stream,indent=2)
        print('clicked only named viewer',a.window,a.x,a.y)
