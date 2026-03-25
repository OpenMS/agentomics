"""Tests for rt_prediction_additive."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestRtPredictionAdditive:
    def test_predict_basic(self):
        from rt_prediction_additive import predict_rt

        result = predict_rt("PEPTIDEK", "krokhin")
        assert result["model"] == "krokhin"
        assert result["length"] == 8
        assert isinstance(result["predicted_rt"], float)

    def test_hydrophobic_peptide_higher_rt(self):
        from rt_prediction_additive import predict_rt

        hydrophobic = predict_rt("LLLLLLLL", "krokhin")
        hydrophilic = predict_rt("KKKKKKKK", "krokhin")
        assert hydrophobic["predicted_rt"] > hydrophilic["predicted_rt"]

    def test_meek_model(self):
        from rt_prediction_additive import predict_rt

        result = predict_rt("PEPTIDEK", "meek")
        assert result["model"] == "meek"
        assert isinstance(result["predicted_rt"], float)

    def test_residue_contributions_sum(self):
        from rt_prediction_additive import predict_rt

        result = predict_rt("PEPTIDEK", "krokhin")
        contrib_sum = sum(c["coefficient"] for c in result["residue_contributions"])
        assert abs(contrib_sum - result["predicted_rt"]) < 0.01

    def test_batch_prediction(self):
        from rt_prediction_additive import predict_batch

        results = predict_batch(["PEPTIDEK", "ANOTHERPEPTIDE"], "krokhin")
        assert len(results) == 2

    def test_different_models_different_results(self):
        from rt_prediction_additive import predict_rt

        krokhin = predict_rt("PEPTIDEK", "krokhin")
        meek = predict_rt("PEPTIDEK", "meek")
        # Different models should give different predictions
        assert krokhin["predicted_rt"] != meek["predicted_rt"]

    def test_residue_contributions_length(self):
        from rt_prediction_additive import predict_rt

        result = predict_rt("ACDEFGHIK", "krokhin")
        assert len(result["residue_contributions"]) == 9
