from paraview.simple import *
from pathlib import Path
import numpy as np
from itertools import product

BLACK = [0., 0., 0.]
WHITE = [1., 1., 1.]

def get_variables_array(max_u, max_c, max_s):
    ALL_CONCENTRATIONS = [
        ["CELLS", "u", max_u],
        ["CELLS", "c_1", max_c],
        ["CELLS", "c_2", max_c],
        ["CELLS", "s_1", max_s],
        ["CELLS", "s_2", max_s],
        ["POINTS", "u_tilde", max_u]
    ]
    return ALL_CONCENTRATIONS

def load_input(input_dir, sim_name, axes_font_size, n_ticks):
    view = GetActiveViewOrCreate("RenderView")
    filename = input_dir + sim_name + ".xdmf"
    print("Loading " + filename)
    reader = Xdmf3ReaderT(registrationName="reader", FileName=filename)
    Show(reader, view=view)
    layout = GetLayout()
    layout.SetSize(800, 800)
    init_config(reader, view, axes_font_size, n_ticks)
    Render()
    return reader, view, layout

def init_config(source, view, axes_font_size=24, n_ticks=4):
    axes = view.AxesGrid
    axes.GridColor = BLACK
    axes.XTitle = axes.YTitle = ""
    axes.XLabelFontFamily = axes.YLabelFontFamily = "Times"
    axes.XLabelFontSize = axes.YLabelFontSize = axes_font_size
    axes.XLabelColor = axes.YLabelColor = BLACK
    axes.AxesToLabel = 3
    axes.FacesToRender = 4
    axes.CullBackface, axes.CullFrontface = 0, 1
    axes.XAxisUseCustomLabels = axes.YAxisUseCustomLabels = 1

    source_bounds = source.GetDataInformation().GetBounds()
    x_min, x_max, y_min, y_max, *_ = source_bounds

    axes.XAxisLabels = np.linspace(x_min, x_max, n_ticks)
    # hard fix for the 0.5-height filter!
    if np.abs(y_max - 0.5) < 0.1:
        _ticks = np.linspace(y_min, 0.5, n_ticks)
        axes.YAxisLabels = np.sort(list(_ticks))
    else:
        axes.YAxisLabels = np.linspace(y_min, y_max, n_ticks)

    axes.UseCustomBounds = 1
    axes.CustomBounds = list(source_bounds)

    view.ResetActiveCameraToNegativeZ()
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Background = WHITE
    view.InteractionMode = "2D"

    y_shift = .075 * abs(y_max - y_min)
    view.ResetCamera(x_min, x_max,
                     y_min - y_shift, y_max - y_shift,
                     0.0, 0.0, True, 0.8)

def postprocess_concentrations(reader, view, rho_P, rho_L):
    # compute additional concentrations c_2, s_2
    c2 = {"Function": "u - c_1",
          "AttributeType": "Cell Data",
          "ResultArrayName": "c_2"}
    _source = Calculator(registrationName="_concentrations",
                         Input=reader, **c2)

    s2 = {"Function": "%f * (1 - u/%f) - s_1" % (rho_L, rho_P),
          "AttributeType": "Cell Data",
          "ResultArrayName": "s_2"};
    source = Calculator(registrationName="concentrations",
                        Input=_source, **s2)
    return source

def plot_contour(source, view, value=1.0, rgb_color=BLACK):
    # plot contour line
    contour_source = Contour(registrationName="u_tilde_contour",
                             Input=source,
                             ContourBy="u_tilde",
                             Isosurfaces=value)
    # config contour appearance
    contour_display = Show(contour_source, view=view)
    contour_display.LineWidth = 3
    contour_display.UseSeparateColorMap = 1
    colormap = GetColorTransferFunction("u_tilde", contour_display)
    colormap.RGBPoints = [0., *rgb_color,
                          1., *rgb_color]
    return contour_source

def plot_velocity(source, view):
    q_array = source.PointData.GetArray("q")
    qx_min, qx_max = q_array.GetRange(0)
    qy_min, qy_max = q_array.GetRange(1)
    q_min = np.sqrt(qx_min**2 + qy_min**2)
    q_max = np.sqrt(qx_max**2 + qy_max**2)
    q_normalizer = 8 * q_max
    quiver_source = Glyph(registrationName="q_quiver",
                          Input=source,
                          GlyphType="Arrow",
                          OrientationArray=["POINTS", "q"],
                          ScaleArray=["POINTS", "q"],
                          ScaleFactor=1/q_normalizer,
                          MaximumNumberOfSamplePoints=150)
    display = GetDisplayProperties(quiver_source, view=view)
    ColorBy(display, ["POINTS", "q"])
    display.SetScalarBarVisibility(view, False)
    colormap = GetColorTransferFunction("q", display)
    colormap.ApplyPreset("Rainbow Desaturated", True)
    return quiver_source

def init_colormaps(display, view, variables,
                  rescale_mode="Never",
                  color_preset="Rainbow Uniform",
                  font_size=24):
    colormaps = colorbars = []
    for i in range(len(variables)):
        var_name = variables[i][1]
        var_max = variables[i][2]

        _map = GetColorTransferFunction(var_name, display,
                                        AutomaticRescaleRangeMode = rescale_mode,
                                        RescaleOnVisibilityChange = 0)
        _map.ApplyPreset(color_preset, True)

        _bar = GetScalarBar(_map, view)
        _bar.Title = _bar.ComponentTitle = ""
        _bar.Orientation = "Horizontal"
        _bar.WindowLocation = "Lower Center"
        _bar.TextPosition = "Ticks left/bottom, annotations right/top"
        _bar.ScalarBarThickness = 16
        _bar.ScalarBarLength = .75
        _bar.LabelColor = BLACK
        _bar.LabelFontFamily = "Times"
        _bar.LabelFontSize = font_size
        _bar.RangeLabelFormat = "%-#6.3g"
        _bar.UseCustomLabels = 1
        _labels = np.linspace(0, var_max, 6)
        _bar.CustomLabels = np.round(_labels)
        _bar.Visibility = 0

        colormaps.append(_map)
        colorbars.append(_bar)
    return colormaps, colorbars

def print_addon(var_name, view,
                suffix, export_path, export_params,
                addon_source, addon_print, addon_name):
    for do_print in addon_print:
        if do_print:
            Show(addon_source, view=view)
            _name = var_name + "_" + addon_name + suffix
        else:
            Hide(addon_source, view=view)
            _name = var_name + suffix
        ExportView(export_path + _name, view=view, **export_params)
    Hide(addon_source, view)

def save_state(out_dir, sim_name):
    state_path = out_dir + sim_name + "/"
    Path(state_path).mkdir(parents=True, exist_ok=True)
    SaveState(state_path + sim_name + "_state.pvsm")

def rescale_colormap(display, array):
    _map = GetColorTransferFunction(array[1], display)
    _map.RescaleTransferFunction(0., array[2])

def export_results(sim_name, input_dir, export_path,
                   export_params, export_times, export_extension,
                   print_bars, print_axes, print_contour, print_quiver,
                   rho_P, rho_L, u_star, max_u, max_c, max_s,
                   axes_font_size, n_ticks):
    reader, view, layout = load_input(input_dir, sim_name,
                                      axes_font_size, n_ticks)
    paraview.simple._DisableFirstRenderCameraReset()

    print("Processing data ...")
    source = postprocess_concentrations(reader, view, rho_P, rho_L)
    contour_args = (plot_contour(source, view, value=u_star),
                    print_contour, "contour")
    quiver_args = (plot_velocity(source, view),
                   print_quiver, "quiver")
    print("Data processed.")

    display = Show(source, view=view)
    ALL_CONCENTRATIONS = get_variables_array(max_u, max_c, max_s)
    init_colormaps(display, view, ALL_CONCENTRATIONS)
    display.RescaleTransferFunctionToDataRange(False, True)
    paraview.simple.HideAll()
    Show(source)

    print("Saving state ...")
    save_state(export_path, sim_name)
    print("Saved state!")

    print("Exporting frames...")
    ext = export_extension
    max_time = max(GetTimeKeeper().TimestepValues)
    for t in export_times:
        print(t)
        view.ViewTime = t
        if t > max_time:
            print("Simulation time exceeded!")
            break
        for j, k in product(print_bars, print_axes):
            suffix = j*"_colorbar" + k*"_axes" + ("_%.0e." % t) + ext
            print(suffix)
            export_args = (suffix, export_path, export_params)
            view.AxesGrid.Visibility = k
            for array in ALL_CONCENTRATIONS:
                ColorBy(display, array[:2])
                display.SetScalarBarVisibility(view, j)
                rescale_colormap(display, array)
                for print_args in [quiver_args, contour_args]:
                    print_addon(array[1], view, *export_args, *print_args)
                display.SetScalarBarVisibility(view, 0)
    print("Frames exported.")
