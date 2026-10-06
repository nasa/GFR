#!MC 1410
#
# Macro to export each frame to a separate png file.
#
$!OpenLayout 'Triangle_Gauss.lay'
#
# Loop over all frames.
$!Loop |NumFrames|
  # Frame names in the layout are 'P0','P1',...,'P8'.
  #
  # Activate Frame.
  $!FrameControl ActivateByNumber
    Frame = |Loop|
  # Export Frame png.
  $!PrintSetup Palette = Color
  $!ExportSetup ExportRegion = CurrentFrame
  # Frame is 1.5 inch square.
  $!ExportSetup ImageWidth = 600
  $!ExportSetup ExportFName = "Triangle_Gauss_|FrameName|.png"
  $!Export 
    ExportRegion = CurrentFrame
$!EndLoop
