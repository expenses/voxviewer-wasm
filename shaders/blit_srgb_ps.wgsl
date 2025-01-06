@binding(1) @group(0) var hdr_0 : texture_2d<f32>;

@binding(3) @group(0) var sampler_linear_clamp_0 : sampler;

@binding(2) @group(0) var tony_mc_mapface_lut_0 : texture_3d<f32>;

struct GlobalParams_std140_0
{
    @align(16) accumulated_frame_index_0 : u32,
};

@binding(0) @group(0) var<uniform> globalParams_0 : GlobalParams_std140_0;
fn tony_mc_mapface_0( stimulus_0 : vec3<f32>) -> vec3<f32>
{
    return (textureSampleLevel((tony_mc_mapface_lut_0), (sampler_linear_clamp_0), (stimulus_0 / (stimulus_0 + vec3<f32>(1.0f)) * vec3<f32>(0.97916668653488159f) + vec3<f32>(0.01041666697710752f)), (0.0f)).xyz);
}

struct pixelOutput_0
{
    @location(0) output_0 : vec4<f32>,
};

struct V2P_0
{
    @builtin(position) Pos_0 : vec4<f32>,
    @location(0) Uv_0 : vec2<f32>,
};

@fragment
fn PSMain( psIn_0 : V2P_0) -> pixelOutput_0
{
    var _S1 : pixelOutput_0 = pixelOutput_0( vec4<f32>(tony_mc_mapface_0((textureSample((hdr_0), (sampler_linear_clamp_0), (psIn_0.Uv_0)).xyz) / vec3<f32>(f32(globalParams_0.accumulated_frame_index_0 + u32(1)))), 1.0f) );
    return _S1;
}

