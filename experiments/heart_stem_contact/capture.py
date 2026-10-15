"""Record actual GUI simulation frames, not a prescribed rotation animation."""
import subprocess
import numpy as np
from PIL import Image


class GuiCapture:
    def __init__(self,folder):
        self.folder=folder; self.frames=0; self.process=None

    def __enter__(self):
        self.log=(self.folder/'ffmpeg.log').open('wb')
        self.process=subprocess.Popen(['ffmpeg','-nostdin','-n','-f','rawvideo','-pixel_format','rgb24',
            '-video_size','960x720','-framerate','30','-i','pipe:0','-an','-c:v','libx264',
            '-preset','fast','-crf','20','-pix_fmt','yuv420p',str(self.folder/'free-contact.mp4')],
            stdin=subprocess.PIPE,stdout=self.log,stderr=self.log)
        return self

    def __call__(self,p,client,t):
        if t+1e-10<self.frames/30: return
        camera=p.getDebugVisualizerCamera(physicsClientId=client)
        frame=p.getCameraImage(960,720,viewMatrix=camera[2],projectionMatrix=camera[3],
            renderer=p.ER_BULLET_HARDWARE_OPENGL,physicsClientId=client)
        rgb=np.asarray(frame[2],dtype=np.uint8).reshape(720,960,4)[:,:,:3].copy()
        if self.frames==0: Image.fromarray(rgb).save(self.folder/'initial.png')
        if self.frames==120: Image.fromarray(rgb).save(self.folder/'at-4s.png')
        self.process.stdin.write(rgb.tobytes()); self.frames+=1

    def __exit__(self,kind,value,traceback):
        self.process.stdin.close(); code=self.process.wait(timeout=30); self.log.close()
        if kind is None and code!=0: raise RuntimeError('video encoder failed')
