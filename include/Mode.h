#pragma once

enum class Mode {
    InEvent = 0, // [OUTDATED]
    OutEvent = 1,
    EndEvent = 2, // [OUTDATED]
    SkipEvent = 3,
    Interaction = 4,
    Init = -1,
};
