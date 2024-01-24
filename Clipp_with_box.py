import vtk
from vtk_read_write import (
   read_poly_data,
   write_poly_data
)

#def read_vtp_file(file_path):
#    reader = vtk.vtkXMLPolyDataReader()
#    reader.SetFileName(file_path)
#    reader.Update()
#
#    return reader.GetOutput()
#
#
#def WritePolyData(input, filename, flag_ascii=0):
#   writer = vtk.vtkXMLPolyDataWriter()
#   writer.SetFileName(filename)
#   writer.SetInputData(input)
#   if flag_ascii:
#    writer.SetDataModeToAscii()
#   writer.Write()


def clip_polyline_with_box(polyData,bounds=(0, 20, -18, 18, -18, 18)):
    # Create a box volume
    box = vtk.vtkBox()
    box.SetBounds(bounds)

    # Clip the Polyline with the box
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(polyData)
    clipper.SetClipFunction(box)
    clipper.SetInsideOut(True)  # Ensure the inside of the box is preserved
    clipper.Update()

    return clipper.GetOutput()


def main():
    # Create and clip the polyline
    polyline = read_poly_data('vtkVilli38Trunc.vtp')
    clipped_polyline = clip_polyline_with_box(polyline)

    # Save the clipped PolyData to a file
    write_poly_data(clipped_polyline, 'vtkVilli38Trunc_clip.vtp',flag_ascii=1)


if __name__ == "__main__":
    main()
