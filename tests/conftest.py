import pytest


@pytest.fixture(scope="function")
def img_data():
    """
    Creating a pytest fixture for a fake Contacts class init.

    :return:
    """
    img_data = []
    return img_data
