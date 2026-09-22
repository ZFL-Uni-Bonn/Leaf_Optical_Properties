from pathlib import Path
import shutil
import re

base_dir = Path("Y:/Wetland Health/Voxelized_Trees")
renamed_dir = base_dir / "renamed"
renamed_dir.mkdir(exist_ok=True)

excluded_folders = {"Voxelized_trees_0.25m"}

file_map = {
    "foliage_voxels.ply": "foliage",
    "wood_voxels.ply": "wood"
}

for voxel_folder in base_dir.glob("Voxelized_trees_*m"):
    if not voxel_folder.is_dir():
        continue

    if voxel_folder.name in excluded_folders:
        print(f"Ausgeschlossen: {voxel_folder.name}")
        continue

    match_res = re.search(r"Voxelized_trees_(.+)m", voxel_folder.name)
    if not match_res:
        continue

    resolution = match_res.group(1)
    target_subfolder = renamed_dir / f"voxelized_trees{resolution}"
    target_subfolder.mkdir(parents=True, exist_ok=True)

    for tree_folder in voxel_folder.iterdir():
        if not tree_folder.is_dir():
            continue

        match_tree = re.match(r"([A-Za-z]+)([0-9.]+)", tree_folder.name)
        if not match_tree:
            print(f"Übersprungen: {tree_folder.name}")
            continue

        species = match_tree.group(1).lower()
        height = match_tree.group(2)

        for old_filename, part_name in file_map.items():
            source_file = tree_folder / old_filename

            if source_file.exists():
                new_filename = f"{species}_{part_name}_h{height}m.ply"
                target_file = target_subfolder / new_filename

                shutil.copy2(source_file, target_file)
                print(f"Kopiert: {source_file} -> {target_file}")
            else:
                print(f"Datei fehlt: {source_file}")