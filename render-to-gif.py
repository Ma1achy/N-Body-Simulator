import imageio.v2 as imageio #imageio v2 is required for gif saving to work properly, do this to get around deprication warning
import os
import datetime

"""
This script will take the frames from a render and save them as a gif

Change the render_dir variable to the directory of the render you want to turn into a gif
Change the starting_frame and end_frame variables to the range of frames you want render as a gif
NOTE: The gif will be saved in the same directory as the frames folder and will be named with the timestamp of when the script was run
      For lots of frames or high resoultion frames, this script will use a lot of memory, so you might want to render a smaller range of frames
"""

def create_gif(frame_dir, gif_path, starting_frame, end_frame):
    """
    save every 4th frame in frame_dir as a gif
    """
    gif_writer = imageio.get_writer(gif_path, duration=playback_speed)

    for frame_number in range(starting_frame, end_frame, 1):  # Increment by 4 to select every 4th frame
        print(f"Creating gif from frame {frame_number} | {((frame_number-starting_frame)/((end_frame-starting_frame)-1)*100):.3f}%", end="\r", flush=True)
        frame_path = os.path.join(frame_dir, f"frame_{frame_number}.png")
        frame = imageio.imread(frame_path)
        gif_writer.append_data(frame)

    print("\nRendering gif...")
    gif_writer.close()

#========================================================================================================

render_dir = r'rendered/2024-01-02 13-16-21 render' #change this to the directory of the render you want to turn into a gif

run_timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S")

frame_dir = os.path.join(render_dir, "frames")
gif_path = os.path.join(render_dir, f"{run_timestamp} render.gif")

playback_speed = 1 #seconds per frame

starting_frame = 0  #frame to start gif at
end_frame = len(os.listdir(frame_dir)) #might not want to use all frames as memory usage is high
#end_frame = 10000

create_gif(frame_dir, gif_path, starting_frame, end_frame)
print(f"\nGIF SAVED TO {gif_path}\n")