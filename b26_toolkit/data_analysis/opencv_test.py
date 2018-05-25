import sys

import cv2
import numpy as np

video = cv2.VideoCapture('Z:\\Lab\\Lev\\lab_measurements\\20171206_sample4_middle_bead\\180114_reset_bead\\cycle_1\\ringdown_1.avi')
video.isOpened()

# Read first frame.
ok, frame = video.read()

# Define an initial bounding box
bbox = (116,80,48,80)

# Uncomment the line below to select a different bounding box
# bbox = cv2.selectROI(frame, False)

tracker = cv2.TrackerMedianFlow_create()

ok = tracker.init(frame, bbox)

while True:
    # Read a new frame
    ok, frame = video.read()
    if not ok:
        break

    # Start timer
    timer = cv2.getTickCount()

    # Update tracker
    ok, bbox = tracker.update(frame)

    # Calculate Frames per second (FPS)
    fps = cv2.getTickFrequency() / (cv2.getTickCount() - timer);

    # Draw bounding box
    if ok:
        # Tracking success
        p1 = (int(bbox[0]), int(bbox[1]))
        p2 = (int(bbox[0] + bbox[2]), int(bbox[1] + bbox[3]))
        cv2.rectangle(frame, p1, p2, (255,0,0), 2, 1)
    else :
        # Tracking failure
        cv2.putText(frame, "Tracking failure detected", (100,80), cv2.FONT_HERSHEY_SIMPLEX, 0.75,(0,0,255),2)

    # # Display tracker type on frame
    # cv2.putText(frame,"Tracker", (100,20), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (50,170,50),2);
    #
    # # Display FPS on frame
    # cv2.putText(frame, "FPS : " + str(int(fps)), (100,50), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (50,170,50), 2);

    print('displaying')
    sys.stdout.flush()

    # Display result
    cv2.imshow("Tracking", frame)

    print('displayed')