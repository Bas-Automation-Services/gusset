/* Gusset, an open sheet metal kernel, as a C library (gusset.dll / libgusset.so).
 *
 * Every function takes a part description as UTF-8 JSON (schema/part.schema.json)
 * and returns a new UTF-8 JSON string, {"ok": ...} or {"error": "..."}, which the
 * caller frees with gusset_free. Lengths in millimetres. */
#ifndef GUSSET_H
#define GUSSET_H

#ifdef __cplusplus
extern "C" {
#endif

/* {"ok": "<flat pattern as DXF R12>"} */
char *gusset_flat_dxf(const char *spec);

/* {"ok": "<flat pattern as SVG>"} */
char *gusset_flat_svg(const char *spec);

/* {"ok": "<folded part as STEP AP214>"} */
char *gusset_step(const char *spec);

/* {"ok": {"volume": mm3, "area": mm2, "center_of_mass": [x, y, z]}} */
char *gusset_mass(const char *spec);

/* Manufacturability findings, default rules:
   {"ok": [{"rule": "...", "feature": "...", "value": mm, "limit": mm, "message": "..."}]} */
char *gusset_check(const char *spec);

/* Frees a string returned above. NULL is ignored. */
void gusset_free(char *text);

#ifdef __cplusplus
}
#endif

#endif
