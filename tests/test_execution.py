import sys
import unittest

import VisHD
from VisHD.config import ExecutionConfig
from VisHD.execution import LocalExecutor, make_executor


def square(value):
    return value * value


class LocalExecutorTests(unittest.TestCase):
    def test_package_import_does_not_load_slurm_dependency(self):
        self.assertIsNotNone(VisHD.load_config)
        self.assertNotIn("submitit", sys.modules)

    def test_submit_returns_result(self):
        job = LocalExecutor().submit(square, 4)
        self.assertEqual(job.wait(), 16)

    def test_map_array_preserves_order(self):
        jobs = LocalExecutor(workers=2).map_array(square, [1, 2, 3])
        self.assertEqual([job.result() for job in jobs], [1, 4, 9])

    def test_factory_defaults_to_local(self):
        executor = make_executor(ExecutionConfig(), "unused")
        self.assertIsInstance(executor, LocalExecutor)


if __name__ == "__main__":
    unittest.main()
