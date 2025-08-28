from sionna.rt import load_scene, Camera

scene = load_scene(
    r"../scenes/office/office_space.xml"
)

my_cam = Camera(position=[9, 35, 0.5], look_at=[0, 0, 3])

scene.render(camera=my_cam)
