from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from VisHD import ConfigError, ExecutionConfig, load_config


class LoadConfigTests(unittest.TestCase):
    def test_local_defaults_and_relative_paths(self):
        with TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            config_path = root / "vishd.toml"
            config_path.write_text(
                """
[paths]
spaceranger_root = "data/spaceranger"
positions_root = "data/positions"
output_root = "results"
""".strip(),
                encoding="utf-8",
            )

            config = load_config(config_path)

            self.assertEqual(config.paths.spaceranger_root, root / "data/spaceranger")
            self.assertEqual(config.paths.positions_root, root / "data/positions")
            self.assertEqual(config.paths.output_root, root / "results")
            self.assertEqual(config.execution, ExecutionConfig())

    def test_slurm_requires_settings(self):
        with TemporaryDirectory() as directory:
            config_path = Path(directory) / "vishd.toml"
            config_path.write_text(
                """
[paths]
spaceranger_root = "data"
output_root = "results"

[execution]
backend = "slurm"
""".strip(),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ConfigError, "execution.slurm is required"):
                load_config(config_path)

    def test_unknown_keys_are_rejected(self):
        with TemporaryDirectory() as directory:
            config_path = Path(directory) / "vishd.toml"
            config_path.write_text(
                """
[paths]
spaceranger_root = "data"
output_root = "results"
typo = "value"
""".strip(),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ConfigError, "unknown key"):
                load_config(config_path)


if __name__ == "__main__":
    unittest.main()
