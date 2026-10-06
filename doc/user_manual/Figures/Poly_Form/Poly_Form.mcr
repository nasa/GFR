#!MC 1410
$!Pick SetMouseMode
  MouseMode = Select
$!Page Name = 'Untitled'
#
# Create a new page.
$!PageControl Create
#
# Load layout file.
$!OpenLayout  "Poly_Form.lay"
#
# Select and export frame 1.
$!FrameControl ActivateByName
  Name = "Frame 001"
$!ExportSetup ExportRegion = CurrentFrame
$!ExportSetup ImageWidth = 1804
$!ExportSetup ExportFName = 'Poly_Form_Modal.png'
$!Export 
  ExportRegion = CurrentFrame
#
# Select and export frame 1.
$!FrameControl ActivateByName
  Name = "Frame 002"
$!ExportSetup ExportFName = 'Poly_Form_Nodal.png'
$!Export 
  ExportRegion = CurrentFrame
