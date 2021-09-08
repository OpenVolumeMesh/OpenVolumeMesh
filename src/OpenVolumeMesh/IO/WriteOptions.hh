#pragma once


struct WriteOptions {
    /// try to determine if a polyhedral mesh can in fact be stored as tet or hex mesh
    bool detect_specialized_mesh = true;
};

