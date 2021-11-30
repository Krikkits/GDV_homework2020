#----------------------------------------------------------
# File __init__.py
#----------------------------------------------------------
 
#    Addon info
bl_info = {
    "name": "RA2 exporter",
    "author": "Manuel Finkh, Benjamin Resch",
    "blender":      (2,6,9),
    "version":      (0,0,1),
    "location": "File > Import-Export",
    "description": "Exporter dumping simply 9 floats for the 3 vertices of each triangle",
    "category": "Import-Export"}
 
# To support reload properly, try to access a package var, 
# if it's there, reload everything
if "bpy" in locals():
    import imp
    if 'ra2_export' in locals():
        imp.reload(ra2_export)
 
import bpy
from bpy.props import *
from bpy_extras.io_utils import ExportHelper, ImportHelper
 
#
#    Export menu
#
 
class EXPORT_OT_ra2(bpy.types.Operator, ExportHelper):
    bl_idname = "io_export_scene.ra2"
    bl_description = 'Export the ra2 format (.ra2)'
    bl_label = "Export ra2"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
 
    # From ExportHelper. Filter filenames.
    filename_ext = ".ra2"
    filter_glob = StringProperty(default="*.ra2", options={'HIDDEN'})
 
    filepath = bpy.props.StringProperty(
        name="File Path", 
        description="File path used for exporting the ra2 file", 
        maxlen= 1024, default= "")
 
    rot90 = bpy.props.BoolProperty(
        name = "Rotate 90 degrees",
        description="Rotate mesh to Y up",
        default = False)
 
    scale = bpy.props.FloatProperty(
        name = "Scale", 
        description="Scale mesh", 
        default = 0.1, min = 0.001, max = 1000.0)
 
    def execute(self, context):
        print("Load", self.properties.filepath)
        from . import ra2_export
        ra2_export.export_ra2(
            self.properties.filepath, 
            context.object, 
            self.rot90, 
            self.scale)
        return {'FINISHED'}
 
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
 
#
#    Registration
#
 
def menu_func_export(self, context):
    self.layout.operator(EXPORT_OT_ra2.bl_idname, text="RA2 file (.ra2)")
 
def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_export.append(menu_func_export)
 
def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_export.remove(menu_func_export)
 
if __name__ == "__main__":
    register()
