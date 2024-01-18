import vtk



def read_vtp_file(file_path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()

    return reader.GetOutput()

def WritePolyData(input,filename,flag_ascii=0):
   writer = vtk.vtkXMLPolyDataWriter()
   writer.SetFileName(filename)
   writer.SetInputData(input)
   if flag_ascii:
    writer.SetDataModeToAscii()
   writer.Write()


def clip_polyline_with_box(polyData):
    # Create a box volume
    box = vtk.vtkBox()
    box.SetBounds(0, 20, -20, 20, -20, 20)

    # Clip the Polyline with the box
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(polyData)
    clipper.SetClipFunction(box)
    clipper.SetInsideOut(True)  # Ensure the inside of the box is preserved
    clipper.Update()

    return clipper.GetOutput()



def main():
    # Create and clip the polyline
    polyline = read_vtp_file('vtkVilli34Trunc.vtp')
    clipped_polyline = clip_polyline_with_box(polyline)

    # Save the clipped PolyData to a file
    WritePolyData(clipped_polyline, 'vtkVilli34Trunc_clip.vtp')

if __name__ == "__main__":
    main()   
