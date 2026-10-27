"""Bounded readers must stop even when the consumer fails with a full queue."""

import threading

import pytest

from arda.annotate.io import prefetched


def test_serial_budget_starts_no_reader():
    caller = threading.get_ident()

    def source():
        assert threading.get_ident() == caller
        yield 1

    with prefetched(source(), enabled=False) as chunks:
        assert list(chunks) == [1]


def test_consumer_failure_stops_the_bounded_reader():
    closed = threading.Event()

    def source():
        try:
            yield from range(10000)
        finally:
            closed.set()

    with pytest.raises(RuntimeError, match="consumer"):
        with prefetched(source()) as chunks:
            assert next(chunks) == 0
            raise RuntimeError("consumer failed")
    assert closed.wait(1)
    assert not any(t.name == "arda-reader" for t in threading.enumerate())


def test_reader_failure_reaches_the_consumer():
    def source():
        yield 1
        raise ValueError("truncated input")

    with pytest.raises(ValueError, match="truncated input"):
        with prefetched(source()) as chunks:
            list(chunks)
