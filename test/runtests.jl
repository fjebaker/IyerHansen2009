using Test
using IyerHansen2009

bs = collect(range(-20, 20, 10))
angles = deflection_angle.(1.0, 0.9, -bs)

# replace NaN for the comparison
replace!(angles, NaN => 0.0)

@test angles ≈ [
    0.22293748909092637,
    0.2964346856086726,
    0.4424230815653054,
    0.874738779448613,
    0.0,
    0.0,
    0.0,
    0.5916923395192697,
    0.3509188629672453,
    0.2511282954376477,
] atol = 1e-9
