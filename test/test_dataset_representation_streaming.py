#!/usr/bin/env python3
"""Regression tests for streamed MolE graph batching."""

from __future__ import annotations

import unittest
from unittest import mock

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data, Dataset, Batch
from rdkit import Chem

from workflow.dataset import dataset_representation as dataset_module


class _FakeGraphDataset(Dataset):
    def __init__(self, size: int, accessed: list[int]) -> None:
        super().__init__()
        self._size = size
        self._accessed = accessed

    def len(self) -> int:
        return self._size

    def get(self, index: int) -> Data:
        self._accessed.append(index)
        return Data(
            x=torch.tensor([[index]], dtype=torch.long),
            edge_index=torch.empty((2, 0), dtype=torch.long),
            edge_attr=torch.empty((0, 2), dtype=torch.long),
            chem_id=f"mol{index}",
        )


class _RecordingModel:
    def __init__(self, accessed: list[int]) -> None:
        self.accessed = accessed
        self.forward_access_snapshots: list[list[int]] = []

    def eval(self) -> "_RecordingModel":
        return self

    def __call__(self, graph_batch: Batch):
        self.forward_access_snapshots.append(list(self.accessed))
        chem_ids = list(graph_batch.chem_id)
        rows = []
        for chem_id in chem_ids:
            index = int(str(chem_id).replace("mol", ""))
            rows.append([float(index), float(index) + 0.5])
        return torch.tensor(rows, dtype=torch.float32), None


def _legacy_representation(dataset: Dataset, model: _RecordingModel, batch_size: int) -> pd.DataFrame:
    graph_list = [dataset[index] for index in range(len(dataset))]
    batch_dataframes = []
    for start in range(0, len(graph_list), batch_size):
        graph_batch = Batch.from_data_list(graph_list[start : start + batch_size])
        representation, _ = model(graph_batch)
        batch_dataframes.append(pd.DataFrame(representation.numpy(), index=graph_batch.chem_id))
    return pd.concat(batch_dataframes)


def _legacy_graph_from_smiles(smiles: str, chem_id: str) -> Data:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    type_idx = []
    chirality_idx = []
    row = []
    col = []
    edge_feat = []

    for atom in mol.GetAtoms():
        type_idx.append(dataset_module.ATOM_INDEX[atom.GetAtomicNum()])
        chirality_idx.append(dataset_module.CHIRALITY_INDEX[atom.GetChiralTag()])

    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        bond_type_idx = dataset_module.BOND_INDEX[bond.GetBondType()]
        bond_dir_idx = dataset_module.BONDDIR_INDEX[bond.GetBondDir()]
        edge_feat.append([bond_type_idx, bond_dir_idx])
        edge_feat.append([bond_type_idx, bond_dir_idx])

    x1 = torch.tensor(type_idx, dtype=torch.long).view(-1, 1)
    x2 = torch.tensor(chirality_idx, dtype=torch.long).view(-1, 1)
    return Data(
        x=torch.cat([x1, x2], dim=-1),
        edge_index=torch.tensor([row, col], dtype=torch.long),
        edge_attr=torch.tensor(np.array(edge_feat), dtype=torch.long),
        chem_id=chem_id,
    )


class BatchRepresentationStreamingTestCase(unittest.TestCase):
    def test_molecule_dataset_matches_legacy_graph_construction_for_representative_smiles(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["CCO", "C[C@H](O)F", "F/C=C/F", "F/C=C\\F", "[He]", "[Na+]"],
                "chem_id": ["achiral", "chiral", "bond_up", "bond_mixed", "helium", "sodium"],
            }
        )

        dataset = dataset_module.MoleculeDataset(smile_df, "smiles", "chem_id")

        for index, row in smile_df.iterrows():
            with self.subTest(smiles=row["smiles"]):
                current = dataset[index]
                legacy = _legacy_graph_from_smiles(row["smiles"], row["chem_id"])
                self.assertEqual(current.chem_id, legacy.chem_id)
                self.assertTrue(torch.equal(current.x, legacy.x))
                self.assertTrue(torch.equal(current.edge_index, legacy.edge_index))
                self.assertTrue(torch.equal(current.edge_attr, legacy.edge_attr))

    def test_molecule_dataset_caches_repeated_raw_smiles_and_preserves_chem_ids(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["C[C@H](O)F", "C[C@H](O)F"],
                "chem_id": ["first", "second"],
            }
        )

        dataset = dataset_module.MoleculeDataset(smile_df, "smiles", "chem_id")

        with mock.patch.object(
            dataset_module.Chem,
            "MolFromSmiles",
            wraps=dataset_module.Chem.MolFromSmiles,
        ) as mol_from_smiles:
            first = dataset[0]
            second = dataset[1]

        self.assertEqual(mol_from_smiles.call_count, 1)
        self.assertEqual(first.chem_id, "first")
        self.assertEqual(second.chem_id, "second")
        self.assertTrue(torch.equal(first.x, second.x))
        self.assertTrue(torch.equal(first.edge_index, second.edge_index))
        self.assertTrue(torch.equal(first.edge_attr, second.edge_attr))

    def test_auto_graph_workers_scale_with_cpu_count_instead_of_fixed_cap(self) -> None:
        with mock.patch.object(dataset_module.os, "cpu_count", return_value=24):
            workers = dataset_module._resolve_graph_worker_count("auto", dataset_size=10_000)

        self.assertEqual(
            workers,
            12,
            "auto graph workers should scale with available CPU capacity instead of staying capped at 4",
        )

    def test_streams_first_graph_batch_before_exhausting_dataset(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN", "CCC", "CCCl", "CCBr"],
                "chem_id": [f"mol{i}" for i in range(5)],
            }
        )
        accessed: list[int] = []
        model = _RecordingModel(accessed)

        with mock.patch.object(
            dataset_module,
            "MoleculeDataset",
            return_value=_FakeGraphDataset(size=5, accessed=accessed),
        ):
            representation = dataset_module.batch_representation(
                smile_df,
                model,
                "smiles",
                "chem_id",
                batch_size=2,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cpu",
            )

        self.assertEqual(model.forward_access_snapshots[0], [0, 1])
        self.assertEqual(list(representation.index), [f"mol{i}" for i in range(5)])
        self.assertEqual(representation.shape, (5, 2))

    def test_iter_batch_representation_yields_embedding_minibatches(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN", "CCC", "CCCl", "CCBr"],
                "chem_id": [f"mol{i}" for i in range(5)],
            }
        )
        accessed: list[int] = []
        model = _RecordingModel(accessed)

        with mock.patch.object(
            dataset_module,
            "MoleculeDataset",
            return_value=_FakeGraphDataset(size=5, accessed=accessed),
        ):
            iterator = dataset_module.iter_batch_representation(
                smile_df,
                model,
                "smiles",
                "chem_id",
                batch_size=2,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cpu",
            )
            first_batch = next(iterator)

        self.assertEqual(accessed, [0, 1])
        self.assertEqual(model.forward_access_snapshots[0], [0, 1])
        self.assertEqual(first_batch["batch_index"], 0)
        self.assertEqual(first_batch["chem_ids"], ["mol0", "mol1"])
        self.assertEqual(list(first_batch["embedding_batch"].index), ["mol0", "mol1"])
        self.assertEqual(first_batch["embedding_batch"].shape, (2, 2))

        remaining = list(iterator)
        self.assertEqual(
            [batch["chem_ids"] for batch in remaining],
            [["mol2", "mol3"], ["mol4"]],
        )
        self.assertEqual(
            [list(batch["embedding_batch"].index) for batch in remaining],
            [["mol2", "mol3"], ["mol4"]],
        )

    def test_iter_batch_representation_yields_batch_profiling(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN", "CCC"],
                "chem_id": [f"mol{i}" for i in range(3)],
            }
        )
        model = _RecordingModel(accessed=[])

        class _ProfiledGraphDataset(_FakeGraphDataset):
            def get(self, index: int) -> Data:
                data = super().get(index)
                data.graph_profile = torch.tensor([[0.1, 0.2, 0.3, 0.4, 1.0]], dtype=torch.float32)
                return data

        with mock.patch.object(
            dataset_module,
            "MoleculeDataset",
            return_value=_ProfiledGraphDataset(size=3, accessed=[]),
        ):
            iterator = dataset_module.iter_batch_representation(
                smile_df,
                model,
                "smiles",
                "chem_id",
                batch_size=2,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cpu",
                enable_profiling=True,
            )
            first_batch = next(iterator)

        profile = first_batch["profiling"]
        self.assertIsInstance(profile, dict)
        self.assertEqual(profile["graph_items"], 2)
        self.assertAlmostEqual(profile["rdkit_parse_seconds"], 0.2, places=6)
        self.assertAlmostEqual(profile["add_hs_seconds"], 0.4, places=6)
        self.assertAlmostEqual(profile["atom_feature_seconds"], 0.6, places=6)
        self.assertAlmostEqual(profile["bond_feature_seconds"], 0.8, places=6)
        self.assertAlmostEqual(profile["graph_total_seconds"], 2.0, places=6)
        self.assertIn("dataloader_iter_seconds", profile)
        self.assertIn("model_forward_seconds", profile)

    def test_streaming_path_preserves_legacy_order_and_values(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN", "CCC", "CCCl", "CCBr"],
                "chem_id": [f"mol{i}" for i in range(5)],
            }
        )

        accessed_streaming: list[int] = []
        accessed_legacy: list[int] = []
        streaming_model = _RecordingModel(accessed_streaming)
        legacy_model = _RecordingModel(accessed_legacy)
        streaming_dataset = _FakeGraphDataset(size=5, accessed=accessed_streaming)
        legacy_dataset = _FakeGraphDataset(size=5, accessed=accessed_legacy)

        with mock.patch.object(dataset_module, "MoleculeDataset", return_value=streaming_dataset):
            streamed = dataset_module.batch_representation(
                smile_df,
                streaming_model,
                "smiles",
                "chem_id",
                batch_size=2,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cpu",
            )

        legacy = _legacy_representation(legacy_dataset, legacy_model, batch_size=2)

        self.assertEqual(list(streamed.index), list(legacy.index))
        self.assertEqual(streamed.shape, legacy.shape)
        np.testing.assert_allclose(streamed.to_numpy(), legacy.to_numpy(), rtol=1e-6, atol=1e-6)

    def test_batch_representation_attaches_profiling_summary_when_enabled(self) -> None:
        smile_df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN", "CCC"],
                "chem_id": [f"mol{i}" for i in range(3)],
            }
        )
        model = _RecordingModel(accessed=[])

        class _ProfiledGraphDataset(_FakeGraphDataset):
            def get(self, index: int) -> Data:
                data = super().get(index)
                data.graph_profile = torch.tensor([[0.1, 0.2, 0.3, 0.4, 1.0]], dtype=torch.float32)
                return data

        with mock.patch.object(
            dataset_module,
            "MoleculeDataset",
            return_value=_ProfiledGraphDataset(size=3, accessed=[]),
        ):
            representation = dataset_module.batch_representation(
                smile_df,
                model,
                "smiles",
                "chem_id",
                batch_size=2,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cpu",
                enable_profiling=True,
            )

        profile = representation.attrs.get("profiling")
        self.assertIsInstance(profile, dict)
        self.assertEqual(profile["graph_items"], 3)
        self.assertAlmostEqual(profile["rdkit_parse_seconds"], 0.3, places=6)
        self.assertAlmostEqual(profile["add_hs_seconds"], 0.6, places=6)
        self.assertAlmostEqual(profile["atom_feature_seconds"], 0.9, places=6)
        self.assertAlmostEqual(profile["bond_feature_seconds"], 1.2, places=6)
        self.assertAlmostEqual(profile["graph_total_seconds"], 3.0, places=6)
        self.assertIn("dataloader_iter_seconds", profile)
        self.assertIn("model_forward_seconds", profile)

    def test_batch_representation_uses_cuda_only_loader_pinning_and_non_blocking_transfer(self) -> None:
        smile_df = pd.DataFrame({"smiles": ["CCO"], "chem_id": ["mol0"]})
        loader_kwargs_log: list[dict] = []
        transfer_log: list[tuple[str, bool]] = []

        class _FakeBatch:
            chem_id = ["mol0"]

            def to(self, device, non_blocking: bool = False):
                transfer_log.append((device, non_blocking))
                return self

        class _FakeLoader:
            def __iter__(self):
                self._yielded = False
                return self

            def __next__(self):
                if self._yielded:
                    raise StopIteration
                self._yielded = True
                return _FakeBatch()

        class _ConstantModel:
            def eval(self):
                return self

            def __call__(self, graph_batch):
                return torch.tensor([[1.0, 2.0]], dtype=torch.float32), None

        def _fake_loader_factory(**kwargs):
            loader_kwargs_log.append(kwargs)
            return _FakeLoader()

        with mock.patch.object(
            dataset_module,
            "MoleculeDataset",
            return_value=_FakeGraphDataset(size=1, accessed=[]),
        ), mock.patch.object(dataset_module, "DataLoader", side_effect=_fake_loader_factory):
            dataset_module.batch_representation(
                smile_df,
                _ConstantModel(),
                "smiles",
                "chem_id",
                batch_size=1,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cpu",
            )
            dataset_module.batch_representation(
                smile_df,
                _ConstantModel(),
                "smiles",
                "chem_id",
                batch_size=1,
                num_graph_workers=0,
                prefetch_batches=2,
                device="cuda:0",
            )

        self.assertFalse(loader_kwargs_log[0].get("pin_memory", False))
        self.assertTrue(loader_kwargs_log[1]["pin_memory"])
        self.assertEqual(transfer_log[0], ("cpu", False))
        self.assertEqual(transfer_log[1], ("cuda:0", True))

    def test_batch_representation_enables_and_restores_deterministic_torch_flags_for_cuda(self) -> None:
        smile_df = pd.DataFrame({"smiles": ["CCO"], "chem_id": ["mol0"]})
        observed_states: list[tuple[bool, bool, bool]] = []

        class _FakeBatch:
            chem_id = ["mol0"]

            def to(self, device, non_blocking: bool = False):
                return self

        class _FakeLoader:
            def __iter__(self):
                self._yielded = False
                return self

            def __next__(self):
                if self._yielded:
                    raise StopIteration
                self._yielded = True
                return _FakeBatch()

        class _DeterminismRecordingModel:
            def eval(self):
                return self

            def __call__(self, graph_batch):
                observed_states.append(
                    (
                        torch.are_deterministic_algorithms_enabled(),
                        torch.backends.cudnn.deterministic,
                        torch.backends.cudnn.benchmark,
                    )
                )
                return torch.tensor([[1.0, 2.0]], dtype=torch.float32), None

        original_algorithms = torch.are_deterministic_algorithms_enabled()
        original_cudnn_deterministic = torch.backends.cudnn.deterministic
        original_cudnn_benchmark = torch.backends.cudnn.benchmark

        torch.use_deterministic_algorithms(False)
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True

        try:
            with mock.patch.object(
                dataset_module,
                "MoleculeDataset",
                return_value=_FakeGraphDataset(size=1, accessed=[]),
            ), mock.patch.object(dataset_module, "DataLoader", return_value=_FakeLoader()):
                dataset_module.batch_representation(
                    smile_df,
                    _DeterminismRecordingModel(),
                    "smiles",
                    "chem_id",
                    batch_size=1,
                    num_graph_workers=0,
                    prefetch_batches=2,
                    device="cuda:0",
                    deterministic_representation=True,
                )
        finally:
            torch.use_deterministic_algorithms(original_algorithms)
            torch.backends.cudnn.deterministic = original_cudnn_deterministic
            torch.backends.cudnn.benchmark = original_cudnn_benchmark

        self.assertEqual(observed_states, [(True, True, False)])
        self.assertEqual(torch.are_deterministic_algorithms_enabled(), original_algorithms)
        self.assertEqual(torch.backends.cudnn.deterministic, original_cudnn_deterministic)
        self.assertEqual(torch.backends.cudnn.benchmark, original_cudnn_benchmark)

    def test_batch_representation_does_not_change_torch_flags_by_default(self) -> None:
        smile_df = pd.DataFrame({"smiles": ["CCO"], "chem_id": ["mol0"]})
        observed_states: list[tuple[bool, bool, bool]] = []

        class _FakeBatch:
            chem_id = ["mol0"]

            def to(self, device, non_blocking: bool = False):
                return self

        class _FakeLoader:
            def __iter__(self):
                self._yielded = False
                return self

            def __next__(self):
                if self._yielded:
                    raise StopIteration
                self._yielded = True
                return _FakeBatch()

        class _DeterminismRecordingModel:
            def eval(self):
                return self

            def __call__(self, graph_batch):
                observed_states.append(
                    (
                        torch.are_deterministic_algorithms_enabled(),
                        torch.backends.cudnn.deterministic,
                        torch.backends.cudnn.benchmark,
                    )
                )
                return torch.tensor([[1.0, 2.0]], dtype=torch.float32), None

        original_algorithms = torch.are_deterministic_algorithms_enabled()
        original_cudnn_deterministic = torch.backends.cudnn.deterministic
        original_cudnn_benchmark = torch.backends.cudnn.benchmark

        torch.use_deterministic_algorithms(False)
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True

        try:
            with mock.patch.object(
                dataset_module,
                "MoleculeDataset",
                return_value=_FakeGraphDataset(size=1, accessed=[]),
            ), mock.patch.object(dataset_module, "DataLoader", return_value=_FakeLoader()):
                dataset_module.batch_representation(
                    smile_df,
                    _DeterminismRecordingModel(),
                    "smiles",
                    "chem_id",
                    batch_size=1,
                    num_graph_workers=0,
                    prefetch_batches=2,
                    device="cuda:0",
                )
        finally:
            torch.use_deterministic_algorithms(original_algorithms)
            torch.backends.cudnn.deterministic = original_cudnn_deterministic
            torch.backends.cudnn.benchmark = original_cudnn_benchmark

        self.assertEqual(observed_states, [(False, False, True)])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
