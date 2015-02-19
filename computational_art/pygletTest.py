import pyglet
import time
import random
import alsaaudio
import audioop

def update_image(dt):
	l,data = inp.read()
	if l:
		index = ((roundint(audioop.rms(data,2), max_volume, num_frames))/(max_volume/num_frames))
		print index
		if index > num_frames-1:
			index = num_frames-1

	img = imgs[index]
	sprite.image = img
	sprite.scale = get_scale(window, img)
	sprite.x = 0
	sprite.y = 0
	window.clear()


window = pyglet.window.Window(fullscreen = True)

@window.event
def on_draw():
	sprite.draw()

def get_scale(window, image):
	if image.width > image.height:
		scale = float(window.width) / image.width
	else:
		scale = float(window.height) / image.height
	return scale

def roundint(val, max_volume, num_frames):
	return int((max_volume/num_frames) * round(float(val)/(max_volume/num_frames)))


if __name__ == '__main__':

	max_volume = 30000
	num_frames = 100

	imgs = [pyglet.image.load('frame-%d.jpg' % i) for i in range(num_frames)]

	inp = alsaaudio.PCM(alsaaudio.PCM_CAPTURE,0)
	inp.setchannels(1)
	inp.setrate(16000)
	inp.setformat(alsaaudio.PCM_FORMAT_S16_LE)
	inp.setperiodsize(160)

	sprites = [pyglet.sprite.Sprite(img) for img in imgs]

	img = random.choice(imgs)
	sprite = pyglet.sprite.Sprite(img)
	sprite.scale = get_scale(window, img) 

	pyglet.clock.schedule_interval(update_image, 1/60.0)

	pyglet.app.run()


