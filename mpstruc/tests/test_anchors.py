from mpstruc.anchors import parse


def test_listing_and_section_anchors_are_read_from_the_page():
    html = ('<a id="id_BETA_BARREL:Outer_Membrane_Autotransporters" ></a> <a id="id_4C00_76_2"></a>'
            '<a id="id_1PRH_1_1"></a> <a id="other"></a> <a id="id_notapdb_1_1"></a>')
    links = parse(html)
    assert links["entries"] == {"4C00": "76_2", "1PRH": "1_1"}
    assert links["sections"] == ["id_BETA_BARREL:Outer_Membrane_Autotransporters"]
    assert links["base"].endswith("/mpstruc/")
