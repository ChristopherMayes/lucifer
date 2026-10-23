# Frames: a run as a series along z

The example runs in one command:

    ../../../production/bin/lucifer lucifer.in

A run writes its beam and its field at element ends. A visualization wants them along z inside the undulators, and it wants to follow a particle from one to the next. This deck turns on the frame series and the field's reductions, and cuts both to four slices of the window ([reading the output](../../doc/reading-output.md)).

The deck runs two segments of the Aramis line with a 1 m comb and a window of 24 wavelengths at 3 wavelengths per slice, which is 8 slices. Each frame carries four of them, indices 3 to 6. A frame of each kind is written at every record, and the record count is in the table below with the run's other measured quantities. A beam frame is a fraction of a megabyte, a field frame about a megabyte, and the statistics file with the reductions a couple of megabytes, so everything the run wrote comes to some tens of megabytes.

The frame index is the statistics record index, so `frames-000006.beam.h5`, `frames-000006.wf.h5` and row 6 of `frames.stats.h5` are one position on the line. Frames land inside the undulators and not only at their ends, which is what riding the comb is for. At a 1 m comb on 8.74 m of line, six of the twelve fall strictly inside an undulator, two on an undulator's end, three in the break between the two, and one at the entry face.

Each file says where it was taken. `sPosition`, `elementName`, `phi0`, `sliceFirst` and `sliceLast` are on every frame. An FEL frame adds `aw`, `ku`, `helical` and `tilt`, and `floorPosition` and `floorAngles` place it in the lab. The averaged mode integrates the quiver away, so a reader that wants the physical orbit rebuilds it from `aw`, `ku` and `s`.

`global%migrate = T` is on so the labels have work to do. Every macroparticle carries a label, and the label follows the particle when migration moves it to a neighbouring slice. Across this series 4591 distinct labels appear in the four slices dumped, against 1024 per slice in any one frame. Particles cross into and out of the range as they slip, and a label that stops appearing has left it. The labels are what make the series trajectories rather than snapshots.

`global%dump_reduced = T` writes the field's intensity projections into the statistics file at every record: `xy_intensity`, `slice_x_intensity`, `slice_y_intensity` and the complex field on axis. They integrate back to `power` exactly and are two hundred times smaller than the field frames, so a picture of the run needs no field frame at all. The frames are here for the views that want the transverse profile of one slice at one position.

The slice range is what keeps the output small. The same deck over the whole window writes twice the bytes, and a run at one wavelength per slice with thousands of slices would write gigabytes a frame. [](../../doc/performance.md) measures a full-size frame at 95 MB and the cost of compression, which does not pay on this data.

The run takes under a second.
