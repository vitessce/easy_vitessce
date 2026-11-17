import pytest
from vitessce.widget import data_server


@pytest.fixture(autouse=True, scope="session")
def stop_vitessce_servers():
    yield
    # Widgets created during tests start uvicorn servers in daemon threads.
    # If they are still running at interpreter shutdown, VitessceWidget.__del__
    # calls thread.join() on a frozen daemon thread, which hangs forever.
    # Stop them while the interpreter is still alive.
    data_server.stop_all()
