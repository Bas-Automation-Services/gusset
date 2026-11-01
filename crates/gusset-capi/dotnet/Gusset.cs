using System.Runtime.InteropServices;
using System.Text;
using System.Text.Json;

namespace GussetKernel;

/// <summary>
/// Gusset's C library (gusset.dll next to the app) from .NET 8: a part description as
/// JSON in (schema/part.schema.json), results out. Millimetres.
/// </summary>
public static class Gusset
{
    const string Dll = "gusset";

    [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)] static extern IntPtr gusset_flat_dxf(byte[] spec);
    [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)] static extern IntPtr gusset_flat_svg(byte[] spec);
    [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)] static extern IntPtr gusset_step(byte[] spec);
    [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)] static extern IntPtr gusset_mass(byte[] spec);
    [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)] static extern IntPtr gusset_check(byte[] spec);
    [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)] static extern void gusset_free(IntPtr text);

    /// <summary>The flat pattern as DXF (R12).</summary>
    public static string FlatDxf(string spec) => Call(gusset_flat_dxf, spec).GetString()!;

    /// <summary>The flat pattern as SVG.</summary>
    public static string FlatSvg(string spec) => Call(gusset_flat_svg, spec).GetString()!;

    /// <summary>The folded part as STEP (AP214).</summary>
    public static string Step(string spec) => Call(gusset_step, spec).GetString()!;

    /// <summary>Volume (mm3), surface area (mm2) and centre of mass (mm) of the folded part.</summary>
    public static (double Volume, double Area, double[] CenterOfMass) Mass(string spec)
    {
        var m = Call(gusset_mass, spec);
        return (m.GetProperty("volume").GetDouble(), m.GetProperty("area").GetDouble(),
            m.GetProperty("center_of_mass").EnumerateArray().Select(x => x.GetDouble()).ToArray());
    }

    /// <summary>Manufacturability findings with the default rules, one message each.</summary>
    public static List<string> Check(string spec) =>
        Call(gusset_check, spec).EnumerateArray().Select(f => f.GetProperty("message").GetString()!).ToList();

    static JsonElement Call(Func<byte[], IntPtr> f, string spec)
    {
        var text = f(Encoding.UTF8.GetBytes(spec + "\0"));
        try
        {
            using var doc = JsonDocument.Parse(Marshal.PtrToStringUTF8(text)!);
            if (doc.RootElement.TryGetProperty("error", out var error))
                throw new GussetException(error.GetString()!);
            return doc.RootElement.GetProperty("ok").Clone();
        }
        finally
        {
            gusset_free(text);
        }
    }
}

/// <summary>The part can't be built: an invalid or (not yet) unsupported description.</summary>
public class GussetException(string message) : Exception(message);
