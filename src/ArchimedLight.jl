module ArchimedLight

import PlantGeom
import PlantGeom: add_ground!
import MultiScaleTreeGraph
import GeometryBasics
import StaticArrays
import LinearAlgebra: norm, dot, cross, det, diag, I
import Serialization
import Dates
import OrderedCollections: OrderedDict
import PlantMeteo
import Adapt
import KernelAbstractions
using KernelAbstractions: @atomic, @index
import Tables
import YAML
import CSV
import JSON
import CodecZlib
import PlyIO

include("types.jl")
include("io.jl")
include("meteo.jl")
include("sky.jl")
include("turtle.jl")
include("interception.jl")
include("scattering.jl")
include("pipeline.jl")
include("voxel/types.jl")
include("voxel/terrain.jl")
include("voxel/io.jl")
include("voxel/terrain_io.jl")
include("voxel/ply_io.jl")
include("voxel/point_input.jl")
include("voxel/traversal_reference.jl")
include("voxel/traversal_dda.jl")
include("voxel/interception.jl")
include("voxel/scattering.jl")
include("voxel/pipeline.jl")
include("voxel/vegetation.jl")
include("voxel/display.jl")
include("voxel/scene.jl")
include("attach.jl")
include("outputs.jl")
include("plantsimengine.jl")
include("visualization.jl")
include("query.jl")
include("voxel/results.jl")

export LightRenderGeometry
export LightNodeMetadata
export LightComponentMetadata
export LightModels
export GroupModel
export TypeModel
export InterceptionModel
export EmitterModel
export OpticalProperties
export LightOptions
export PixelHitStackPolicy
export AutoPixelHitStack
export SmallPixelHitStack
export VectorPixelHitStack
export RasterCPUBackend, RasterGPUBackendConfig, RasterGPUBackend
export RaycastScatteringBackend, RasterGPUScatteringBackend
export SkyState
export LightBudget
export LightStepResult
export ValidationReport
export SceneSummary
export MeteoSummary
export LightSimulation
export VoxelGrid
export VoxelScene
export add_ground!
export VoxelPointInput, VoxelPointComponent
export VoxelMaterial, VoxelVegetation, VoxelVegetationStepResult, VoxelComponentLight
export AbstractVoxelTerrain
export NoVoxelTerrain
export TerrainHit
export SoilOpticalProperties
export PlanarTerrain
export HeightFieldTerrain
export TriangulatedTerrain
export TerrainCoordinateTransform
export AAIGrid
export GroundDistanceAlignment
export VoxelTerrainSceneInput
export VoxelTerrainScene
export VoxelRaySegment
export VoxelRayPath
export VoxelCPUBackend
export VoxelOpticalProperties
export VoxelGroundOptics
export VoxelDirectionResponse
export VoxelResponseCache
export VoxelFirstOrderResult
export VoxelScatteringQuadrature
export VoxelScatteringTransportCache
export VoxelTransportResult
export VoxelBandScatteringResult
export VoxelScatteringResult
export VoxelLightStepResult

export read_scene
export write_scene
export read_models
export prepare_models
export models_for
export translucent
export virtual_sensor
export emitter
export read_options
export read_simulation
export read_meteo
export prepare_meteo
export interpolate_radiation
export cache_summary
export update_scene!
export update_models!
export update_options!
export check_scene
export check_models
export check_meteo
export check_simulation
export summarize_scene
export summarize_meteo
export run_light
export read_voxel_grid
export read_voxel_point_input, read_ply_terrain, voxel_component_values
export read_voxel_column
export write_voxel_grid
export read_aai_grid
export height_field_terrain
export validate_ground_distance
export terrain_from_ground_distance
export read_voxel_terrain_scene
export intersect_terrain
export terrain_bounds
export terrain_patch_count
export terrain_patch_point
export terrain_patch_normal
export terrain_patch_material
export terrain_elevation_at
export terrain_lambertian_weights
export soil_optics
export soil_absorbed_energy
export soil_reflected_energy
export trace_voxel_ray
export prepare_voxel_responses
export compute_voxel_first_order
export prepare_voxel_scattering_quadrature
export prepare_voxel_scattering_transport
export apply_voxel_scattering_transport
export compute_voxel_scattering
export voxel_generic_green_leaf_optics
export voxel_single_scattering_albedo
export voxel_absorptance
export VOXEL_PAR_WAVELENGTH_NM
export VOXEL_NIR_WAVELENGTH_NM
export run_voxel_light_step
export run_voxel_light_series
export attach_node_values!
export attach_light_step!
export attach_light_series!
export write_component_values
export CompiledComponentAggregation
export component_values
export component_values!
export ArchimedLightModel
export archimed_light_outputs
export light_render_geometry
export tile_light_geometry
export light_node_ids
export light_metric_values
export light_face_values
export light_vertex_values
export lightplot
export lightplot!

end
