/*
    Visualization tool for xyz files
    Uses raylib (version 5.5.0)
    (written by claude)
    coordinates are in reduced units (sigma_Ar = 1); each sphere is drawn
    with radius sigma_i / 2 from elements.cpp, so particles at distance
    sigma_ij just touch

    Usage:  bin/viewer [file1.xyz file2.xyz ...] [font.ttf]
            (default file: input/test.xyz)
            Several files are played as an animation, in the order given:
                bin/viewer output/initialconfig.xyz $(ls -v output/step*.xyz)
            (ls -v sorts step2 before step10). An argument ending in .ttf or
            .otf is used as the font.

    Mouse:  left drag   rotate around the particles (or move the frame slider)
            right drag  pan
            wheel       zoom
    Keys:   V           toggle velocity arrows
            + / -       longer / shorter velocity arrows
            B           toggle simulation box
            R           reset view
            Space       play / pause the animation
            Left/Right  previous / next frame (pauses)
            Up/Down     faster / slower playback
            Home or 0   back to frame 0
    Frames are numbered from 0 (= first file, e.g. initialconfig.xyz).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "raylib.h"
#include "xyzio.h"
#include "elements.h"
#include "paths.h"

#define DEFAULT_FONT "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
#define FONT_SIZE    22     // size the text is drawn at (pixels)
#define PLAY_FPS     10.0f  // default playback speed (trajectory frames per second)

// ---------------------------------------------------------------------
// Colors per label. Known elements get fixed colors (roughly the usual
// CPK/Jmol scheme), any other label gets the next color from the palette
// in order of first appearance.
// ---------------------------------------------------------------------
typedef struct { const char *label; Color color; } LabelColor;

static const LabelColor KNOWN[] = {
    { "H",  { 230, 230, 230, 255 } }, { "C",  {  80,  80,  80, 255 } },
    { "N",  {  48,  80, 248, 255 } }, { "O",  { 255,  13,  13, 255 } },
    { "He", { 217, 255, 255, 255 } }, { "Ne", { 179, 227, 245, 255 } },
    { "Ar", { 128, 209, 227, 255 } }, { "Kr", {  92, 184, 209, 255 } },
    { "Xe", {  66, 158, 176, 255 } },
};

static const Color PALETTE[] = {
    { 230, 159,   0, 255 }, {  86, 180, 233, 255 }, {   0, 158, 115, 255 },
    { 213,  94,   0, 255 }, { 204, 121, 167, 255 }, {   0, 114, 178, 255 },
    { 240, 228,  66, 255 }, { 120, 120, 120, 255 },
};

#define MAX_TYPES 32

typedef struct {
    int  n;
    char label[MAX_TYPES][LABEL_LEN];
    Color color[MAX_TYPES];
} TypeTable;

// Returns the index of `label` in the table, adding it if it is new.
static int type_index(TypeTable *t, const char *label) {
    for (int k = 0; k < t->n; ++k)
        if (strcmp(t->label[k], label) == 0) return k;
    if (t->n == MAX_TYPES) return MAX_TYPES - 1;   // too many types: share the last color

    int k = t->n++;
    snprintf(t->label[k], LABEL_LEN, "%s", label);
    static int next = 0;
    t->color[k] = PALETTE[next++ % (sizeof(PALETTE) / sizeof(PALETTE[0]))];
    for (size_t j = 0; j < sizeof(KNOWN) / sizeof(KNOWN[0]); ++j)
        if (strcmp(KNOWN[j].label, label) == 0) { t->color[k] = KNOWN[j].color; --next; break; }
    return k;
}

// raylib works in float, the frames in double: convert explicitly
static Vector3 v3(double x, double y, double z) {
    return { (float)x, (float)y, (float)z };
}
static Vector3 v3(Vec3 r) { return v3(r.x, r.y, r.z); }

// ---------------------------------------------------------------------
// 3D arrow from the surface of the sphere at p (radius r) along v.
// Length is scale*|v|: a cylinder shaft plus a cone head.
// ---------------------------------------------------------------------
static void draw_arrow(Vector3 p, float r, Vector3 v, float scale, Color col) {
    float vl = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
    float len = scale * vl;
    if (len < 1e-3f) return;                       // (near) zero velocity: no arrow
    Vector3 u = { v.x / vl, v.y / vl, v.z / vl };  // unit direction

    const float shaft_r = 0.035f, head_r = 0.10f;
    float head_len = fminf(0.25f, 0.6f * len);     // short arrows are mostly head

    Vector3 a = { p.x + r * u.x, p.y + r * u.y, p.z + r * u.z };           // start
    Vector3 b = { a.x + (len - head_len) * u.x, a.y + (len - head_len) * u.y,
                  a.z + (len - head_len) * u.z };                          // head base
    Vector3 tip = { a.x + len * u.x, a.y + len * u.y, a.z + len * u.z };

    if (len > head_len) DrawCylinderEx(a, b, shaft_r, shaft_r, 12, col);
    DrawCylinderEx(b, tip, head_r, 0.0f, 16, col);
}

// Orbit camera: the camera sits on a sphere of radius `dist` around `target`,
// at angles yaw (around the y axis) and pitch (up/down).
typedef struct {
    Vector3 target;
    float dist, yaw, pitch;
} Orbit;

static void orbit_apply(const Orbit *o, Camera3D *cam) {
    cam->target = o->target;
    cam->position = Vector3{
        o->target.x + o->dist * cosf(o->pitch) * sinf(o->yaw),
        o->target.y + o->dist * sinf(o->pitch),
        o->target.z + o->dist * cosf(o->pitch) * cosf(o->yaw)
    };
}

static void orbit_update(Orbit *o, const Camera3D *cam) {
    Vector2 d = GetMouseDelta();

    // rotate
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        o->yaw   -= 0.005f * d.x;
        o->pitch += 0.005f * d.y;
        const float lim = 1.55f;    // just under 90 deg, avoids flipping over the pole
        if (o->pitch >  lim) o->pitch =  lim;
        if (o->pitch < -lim) o->pitch = -lim;
    }

    // zoom (multiplicative, so it feels the same near and far)
    float wheel = GetMouseWheelMove();
    if (wheel != 0.0f) {
        o->dist *= powf(0.9f, wheel);
        if (o->dist < 0.5f) o->dist = 0.5f;
    }

    // pan: move the target in the camera's screen plane
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
        Vector3 fw = { cam->target.x - cam->position.x,
                       cam->target.y - cam->position.y,
                       cam->target.z - cam->position.z };
        float fl = sqrtf(fw.x*fw.x + fw.y*fw.y + fw.z*fw.z);
        fw.x /= fl; fw.y /= fl; fw.z /= fl;

        // right = fw x (0,1,0),  up = right x fw
        Vector3 rt = { -fw.z, 0.0f, fw.x };
        float rl = sqrtf(rt.x*rt.x + rt.z*rt.z);
        rt.x /= rl; rt.z /= rl;
        Vector3 up = { rt.y*fw.z - rt.z*fw.y,
                       rt.z*fw.x - rt.x*fw.z,
                       rt.x*fw.y - rt.y*fw.x };

        float s = 0.0015f * o->dist;     // pan faster when zoomed out
        o->target.x += s * (-d.x * rt.x + d.y * up.x);
        o->target.y += s * (-d.x * rt.y + d.y * up.y);
        o->target.z += s * (-d.x * rt.z + d.y * up.z);
    }
}

int main(int argc, char **argv) {
    // split the arguments: .ttf/.otf -> font, everything else -> trajectory frames
    const char *font_path = DEFAULT_FONT;
    const char **paths = static_cast<const char**>(malloc((argc > 1 ? argc : 2) * sizeof *paths));
    int nframes = 0;
    for (int a = 1; a < argc; ++a) {
        if (IsFileExtension(argv[a], ".ttf;.otf")) font_path = argv[a];
        else paths[nframes++] = argv[a];
    }
    if (nframes == 0) paths[nframes++] = INPUT_DIR "test.xyz";

    // load every frame up front; all frames must describe the same particles
    XYZFrame **frames = static_cast<XYZFrame**>(calloc(nframes, sizeof *frames));
    for (int k = 0; k < nframes; ++k) {
        frames[k] = xyz_read(paths[k]);       // prints the reason on failure
        int ok = frames[k] != NULL;
        if (ok && frames[k]->n() != frames[0]->n()) {
            fprintf(stderr, "viewer: %s has %d particles, %s has %d\n",
                    paths[k], frames[k]->n(), paths[0], frames[0]->n());
            ok = 0;
        }
        for (int i = 0; ok && i < frames[k]->n(); ++i)
            if (frames[k]->type(i) != frames[0]->type(i)) {
                fprintf(stderr, "viewer: particle %d is '%s' in %s but '%s' in %s\n",
                        i, frames[k]->label(i), paths[k], frames[0]->label(i), paths[0]);
                ok = 0;
            }
        if (!ok) {
            for (int j = 0; j <= k; ++j) xyz_free(frames[j]);
            free(frames); free(paths);
            return 1;
        }
    }
    XYZFrame *f = frames[0];                  // frame currently shown

    // camera centered on the first frame's box
    Vector3 c = v3(f->box().x / 2, f->box().y / 2, f->box().z / 2);
    float Lmax = (float)fmax(f->box().x, fmax(f->box().y, f->box().z));

    SetConfigFlags(FLAG_MSAA_4X_HINT);    // smoother sphere edges
    InitWindow(1000, 800, nframes == 1 ? GetFileName(paths[0])            // title: file name without folders
                                       : TextFormat("%s ... %s (%d frames)", GetFileName(paths[0]),
                                                    GetFileName(paths[nframes - 1]), nframes));
    SetTargetFPS(60);

    // Load the font at twice the drawn size and filter it, so it stays sharp.
    // Falls back to raylib's built-in pixel font if the file can't be loaded.
    Font font = GetFontDefault();
    if (FileExists(font_path)) {
        Font tmp = LoadFontEx(font_path, 2 * FONT_SIZE, NULL, 0);
        if (IsFontValid(tmp)) {
            font = tmp;
            SetTextureFilter(font.texture, TEXTURE_FILTER_BILINEAR);
        }
    } else {
        fprintf(stderr, "viewer: font '%s' not found, using default font\n", font_path);
    }

    const Orbit home = { .target = c, .dist = 2.5f * Lmax, .yaw = 0.785f, .pitch = 0.615f };
    Orbit orbit = home;
    Camera3D cam = { .position = {}, .target = {}, .up = { 0, 1, 0 }, .fovy = 45, .projection = CAMERA_PERSPECTIVE };
    orbit_apply(&orbit, &cam);

    float vscale = 0.5f;          // arrow length = vscale * |v|, change with + / -
    int show_v = 1;
    int show_box = 1;

    // animation state
    int   cur = 0;                  // index of the frame shown
    int   playing = nframes > 1;
    float fps = PLAY_FPS;
    float acc = 0.0f;               // time collected towards the next frame
    int   dragging = 0;             // 1 while the frame slider is held with the mouse

    // assign each particle a type (and thus a color) once, up front
    TypeTable types{};
    std::vector<int> type_of(f->n());
    for (int i = 0; i < f->n(); ++i) type_of[i] = type_index(&types, f->label(i));

    while (!WindowShouldClose()) {
        // frame slider along the bottom, above the help lines
        Rectangle track = { 20, (float)(GetScreenHeight() - 92), (float)(GetScreenWidth() - 40), 6 };
        Rectangle hit   = { track.x - 10, track.y - 12, track.width + 20, track.height + 24 };
        if (nframes > 1) {
            Vector2 m = GetMousePosition();
            if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(m, hit)) {
                dragging = 1;
                playing = 0;
            }
            if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) dragging = 0;
            if (dragging) {
                float t = fminf(fmaxf((m.x - track.x) / track.width, 0.0f), 1.0f);
                cur = (int)lroundf(t * (nframes - 1));
            }
        }

        // the left mouse button either drags the slider or rotates the camera, not both
        if (!dragging) orbit_update(&orbit, &cam);
        if (IsKeyPressed(KEY_R)) orbit = home;
        if (IsKeyPressed(KEY_V)) show_v = !show_v;
        // +/- by typed character, so it works on any keyboard layout
        for (int ch = GetCharPressed(); ch != 0; ch = GetCharPressed()) {
            if (ch == '+') vscale = fminf(vscale * 1.5f, 10.0f);
            if (ch == '-') vscale = fmaxf(vscale / 1.5f, 0.02f);
            if (ch == '0') { cur = 0; acc = 0.0f; }
        }
        if (IsKeyPressed(KEY_HOME)) { cur = 0; acc = 0.0f; }
        if (IsKeyPressed(KEY_B)) show_box = !show_box;
        orbit_apply(&orbit, &cam);

        // animation: keys, then advance by elapsed time (independent of the 60 FPS redraw)
        if (nframes > 1) {
            if (IsKeyPressed(KEY_SPACE)) { playing = !playing; acc = 0.0f; }
            if (IsKeyPressed(KEY_RIGHT)) { playing = 0; cur = (cur + 1) % nframes; }
            if (IsKeyPressed(KEY_LEFT))  { playing = 0; cur = (cur - 1 + nframes) % nframes; }
            if (IsKeyPressed(KEY_UP))    fps = fminf(fps * 1.5f, 60.0f);
            if (IsKeyPressed(KEY_DOWN))  fps = fmaxf(fps / 1.5f, 0.5f);
            if (playing) {
                acc += GetFrameTime();
                while (acc >= 1.0f / fps) { acc -= 1.0f / fps; cur = (cur + 1) % nframes; }
            }
        }
        f = frames[cur];
        Vector3 L = v3(f->box());
        Vector3 bc = { L.x / 2, L.y / 2, L.z / 2 };   // box center of this frame

        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode3D(cam);
            DrawGrid(20, 1.0f);
            if (show_box) DrawCubeWires(bc, L.x, L.y, L.z, DARKGRAY);
            for (int i = 0; i < f->n(); ++i) {
                Vector3 p = v3(f->pos(i));
                Color col = types.color[type_of[i]];
                float radius = (float)element_radius(f->label(i));   // sigma_i / 2
                DrawSphere(p, radius, col);
                DrawSphereWires(p, radius, 8, 8, ColorBrightness(col, -0.4f));
                if (show_v) {
                    Vector3 v = v3(f->vel(i));
                    draw_arrow(p, radius, v, vscale, Color{ 200, 30, 30, 255 });
                }
            }
        EndMode3D();
        DrawTextEx(font, f->comment.c_str(), Vector2{ 10, 10 }, FONT_SIZE, 1, DARKGRAY);
        DrawTextEx(font, TextFormat("N = %d   box: %g x %g x %g (%s)", f->n(), L.x, L.y, L.z, boundary_name(f->boundary())), Vector2{ 10, 38 }, FONT_SIZE, 1, DARKGRAY);

        // legend: one colored dot + label per particle type
        for (int k = 0; k < types.n; ++k) {
            float y = 76 + k * (FONT_SIZE + 6);
            DrawCircle(20, (int)y + FONT_SIZE / 2, 8, types.color[k]);
            DrawCircleLines(20, (int)y + FONT_SIZE / 2, 8, ColorBrightness(types.color[k], -0.4f));
            DrawTextEx(font, types.label[k], Vector2{ 36, y }, FONT_SIZE, 1, DARKGRAY);
        }
        if (nframes > 1) {
            const char *status = TextFormat("frame %d / %d   %s   %.1f frames/s", cur, nframes - 1,
                                            playing ? "playing" : "paused", fps);
            float w = MeasureTextEx(font, status, FONT_SIZE, 1).x;
            DrawTextEx(font, status, Vector2{ GetScreenWidth() - w - 10, 10 }, FONT_SIZE, 1, DARKGRAY);
            DrawTextEx(font, "[Space] play/pause   [Left/Right] step   [Up/Down] speed   [Home/0] frame 0",
                       Vector2{ 10, (float)(GetScreenHeight() - 54) }, 18, 1, GRAY);

            // slider: track, played part, a tick per frame (if they fit), knob
            float kx = track.x + track.width * cur / (float)(nframes - 1);
            float cy = track.y + track.height / 2;
            DrawRectangleRounded(track, 1.0f, 6, LIGHTGRAY);
            DrawRectangleRounded(Rectangle{ track.x, track.y, kx - track.x, track.height }, 1.0f, 6, SKYBLUE);
            if (track.width / (nframes - 1) >= 6)
                for (int k = 0; k < nframes; ++k) {
                    float tx = track.x + track.width * k / (float)(nframes - 1);
                    DrawLineEx(Vector2{ tx, track.y + track.height + 3 },
                               Vector2{ tx, track.y + track.height + 8 }, 1.0f, GRAY);
                }
            Vector2 mouse = GetMousePosition();
            int hover = dragging || CheckCollisionPointRec(mouse, hit);
            DrawCircleV(Vector2{ kx, cy }, hover ? 10 : 8, hover ? DARKBLUE : BLUE);
            DrawCircleLinesV(Vector2{ kx, cy }, hover ? 10 : 8, RAYWHITE);
            const char *lbl = TextFormat("%d", cur);
            float lw = MeasureTextEx(font, lbl, 16, 1).x;
            DrawTextEx(font, lbl, Vector2{ kx - lw / 2, track.y - 30 }, 16, 1, DARKGRAY);
        }
        DrawTextEx(font, "left drag: rotate   right drag: pan   wheel: zoom   [R] reset   [V] velocities   [+/-] arrow length   [B] box",
                   Vector2{ 10, (float)(GetScreenHeight() - 30) }, 18, 1, GRAY);
        EndDrawing();
    }

    UnloadFont(font);     // no-op for the default font
    CloseWindow();
    for (int k = 0; k < nframes; ++k) xyz_free(frames[k]);
    free(frames); free(paths);
    return 0;
}
